// Copyright (c) 2018, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

#include "estimators/two_view_geometry.h"

#include <unordered_set>

#include "base/camera.h"
#include "base/essential_matrix.h"
#include "base/homography_matrix.h"
#include "base/pose.h"
#include "base/projection.h"
#include "base/triangulation.h"
#include "estimators/essential_matrix.h"
#include "estimators/fundamental_matrix.h"
#include "estimators/homography_matrix.h"
#include "estimators/translation_transform.h"
#include "optim/loransac.h"
#include "optim/ransac.h"
#include "util/random.h"

// IS: GCRansac additions
#include "GCRANSAC.h"
#include "essential_estimator.h"
#include "fundamental_estimator.h"
#include "gcransac_utils.h"
#include "grid_neighborhood_graph.h"
#include "homography_estimator.h"
#include "model.h"
#include "progressive_napsac_sampler.h"
#include "uniform_sampler.h"
#include <opencv2/core/core.hpp>

namespace colmap {
namespace {

FeatureMatches ExtractInlierMatches(const FeatureMatches& matches,
                                    const size_t num_inliers,
                                    const std::vector<char>& inlier_mask) {
  FeatureMatches inlier_matches(num_inliers);
  size_t j = 0;
  for (size_t i = 0; i < matches.size(); ++i) {
    if (inlier_mask[i]) {
      inlier_matches[j] = matches[i];
      j += 1;
    }
  }
  return inlier_matches;
}

FeatureMatches ExtractOutlierMatches(const FeatureMatches& matches,
                                     const FeatureMatches& inlier_matches) {
  CHECK_GE(matches.size(), inlier_matches.size());

  std::unordered_set<std::pair<point2D_t, point2D_t>> inlier_matches_set;
  inlier_matches_set.reserve(inlier_matches.size());
  for (const auto& match : inlier_matches) {
    inlier_matches_set.emplace(match.point2D_idx1, match.point2D_idx2);
  }

  FeatureMatches outlier_matches;
  outlier_matches.reserve(matches.size() - inlier_matches.size());

  for (const auto& match : matches) {
    if (inlier_matches_set.count(
            std::make_pair(match.point2D_idx1, match.point2D_idx2)) == 0) {
      outlier_matches.push_back(match);
    }
  }

  return outlier_matches;
}

inline bool IsImagePointInBoundingBox(const Eigen::Vector2d& point,
                                      const double minx, const double maxx,
                                      const double miny, const double maxy) {
  return point.x() >= minx && point.x() <= maxx && point.y() >= miny &&
         point.y() <= maxy;
}

}  // namespace

void TwoViewGeometry::Invert() {
  F.transposeInPlace();
  E.transposeInPlace();
  H = H.inverse().eval();

  const Eigen::Vector4d orig_qvec = qvec;
  const Eigen::Vector3d orig_tvec = tvec;
  InvertPose(orig_qvec, orig_tvec, &qvec, &tvec);

  for (auto& match : inlier_matches) {
    std::swap(match.point2D_idx1, match.point2D_idx2);
  }
}

void TwoViewGeometry::Estimate(const Camera& camera1,
                               const std::vector<Eigen::Vector2d>& points1,
                               const Camera& camera2,
                               const std::vector<Eigen::Vector2d>& points2,
                               const FeatureMatches& matches,
                               const Options& options) {
  if (camera1.HasPriorFocalLength() && camera2.HasPriorFocalLength()) {
    EstimateCalibrated(camera1, points1, camera2, points2, matches, options);
  } else {
    EstimateUncalibrated(camera1, points1, camera2, points2, matches, options);
  }
}

void TwoViewGeometry::EstimateMultiple(
    const Camera& camera1, const std::vector<Eigen::Vector2d>& points1,
    const Camera& camera2, const std::vector<Eigen::Vector2d>& points2,
    const FeatureMatches& matches, const Options& options) {
  FeatureMatches remaining_matches = matches;
  std::vector<TwoViewGeometry> two_view_geometries;
  while (true) {
    TwoViewGeometry two_view_geometry;
    two_view_geometry.Estimate(camera1, points1, camera2, points2,
                               remaining_matches, options);
    if (two_view_geometry.config == ConfigurationType::DEGENERATE) {
      break;
    }

    if (options.multiple_ignore_watermark) {
      if (two_view_geometry.config != ConfigurationType::WATERMARK) {
        two_view_geometries.push_back(two_view_geometry);
      }
    } else {
      two_view_geometries.push_back(two_view_geometry);
    }

    remaining_matches = ExtractOutlierMatches(remaining_matches,
                                              two_view_geometry.inlier_matches);
  }

  if (two_view_geometries.empty()) {
    config = ConfigurationType::DEGENERATE;
  } else if (two_view_geometries.size() == 1) {
    *this = two_view_geometries[0];
  } else {
    config = ConfigurationType::MULTIPLE;

    for (const auto& two_view_geometry : two_view_geometries) {
      inlier_matches.insert(inlier_matches.end(),
                            two_view_geometry.inlier_matches.begin(),
                            two_view_geometry.inlier_matches.end());
    }
  }
}

bool TwoViewGeometry::EstimateRelativePose(
    const Camera& camera1, const std::vector<Eigen::Vector2d>& points1,
    const Camera& camera2, const std::vector<Eigen::Vector2d>& points2) {
  // We need a valid epopolar geometry to estimate the relative pose.
  if (config != CALIBRATED && config != UNCALIBRATED && config != PLANAR &&
      config != PANORAMIC && config != PLANAR_OR_PANORAMIC) {
    return false;
  }

  // Extract normalized inlier points.
  std::vector<Eigen::Vector2d> inlier_points1_normalized;
  inlier_points1_normalized.reserve(inlier_matches.size());
  std::vector<Eigen::Vector2d> inlier_points2_normalized;
  inlier_points2_normalized.reserve(inlier_matches.size());
  for (const auto& match : inlier_matches) {
    const point2D_t idx1 = match.point2D_idx1;
    const point2D_t idx2 = match.point2D_idx2;
    inlier_points1_normalized.push_back(camera1.ImageToWorld(points1[idx1]));
    inlier_points2_normalized.push_back(camera2.ImageToWorld(points2[idx2]));
  }

  Eigen::Matrix3d R;
  std::vector<Eigen::Vector3d> points3D;

  if (config == CALIBRATED || config == UNCALIBRATED) {
    // Try to recover relative pose for calibrated and uncalibrated
    // configurations. In the uncalibrated case, this most likely leads to a
    // ill-defined reconstruction, but sometimes it succeeds anyways after e.g.
    // subsequent bundle-adjustment etc.
    PoseFromEssentialMatrix(E, inlier_points1_normalized,
                            inlier_points2_normalized, &R, &tvec, &points3D);
  } else if (config == PLANAR || config == PANORAMIC ||
             config == PLANAR_OR_PANORAMIC) {
    Eigen::Vector3d n;
    PoseFromHomographyMatrix(
        H, camera1.CalibrationMatrix(), camera2.CalibrationMatrix(),
        inlier_points1_normalized, inlier_points2_normalized, &R, &tvec, &n,
        &points3D);
  } else {
    return false;
  }

  qvec = RotationMatrixToQuaternion(R);

  if (points3D.empty()) {
    tri_angle = 0;
  } else {
    tri_angle = Median(CalculateTriangulationAngles(
        Eigen::Vector3d::Zero(), -R.transpose() * tvec, points3D));
  }

  if (config == PLANAR_OR_PANORAMIC) {
    if (tvec.norm() == 0) {
      config = PANORAMIC;
      tri_angle = 0;
    } else {
      config = PLANAR;
    }
  }

  return true;
}

void TwoViewGeometry::EstimateCalibrated(
    const Camera& camera1, const std::vector<Eigen::Vector2d>& points1,
    const Camera& camera2, const std::vector<Eigen::Vector2d>& points2,
    const FeatureMatches& matches, const Options& options) {
  options.Check();

  if (matches.size() < options.min_num_inliers) {
    config = ConfigurationType::DEGENERATE;
    return;
  }

  // Extract corresponding points.
  std::vector<Eigen::Vector2d> matched_points1(matches.size());
  std::vector<Eigen::Vector2d> matched_points2(matches.size());
  std::vector<Eigen::Vector2d> matched_points1_normalized(matches.size());
  std::vector<Eigen::Vector2d> matched_points2_normalized(matches.size());
  for (size_t i = 0; i < matches.size(); ++i) {
    const point2D_t idx1 = matches[i].point2D_idx1;
    const point2D_t idx2 = matches[i].point2D_idx2;
    matched_points1[i] = points1[idx1];
    matched_points2[i] = points2[idx2];
    matched_points1_normalized[i] = camera1.ImageToWorld(points1[idx1]);
    matched_points2_normalized[i] = camera2.ImageToWorld(points2[idx2]);
  }

  // Estimate epipolar models.

  auto E_ransac_options = options.ransac_options;
  E_ransac_options.max_error =
      (camera1.ImageToWorldThreshold(options.ransac_options.max_error) +
       camera2.ImageToWorldThreshold(options.ransac_options.max_error)) /
      2;

  LORANSAC<EssentialMatrixFivePointEstimator, EssentialMatrixFivePointEstimator>
      E_ransac(E_ransac_options);
  const auto E_report =
      E_ransac.Estimate(matched_points1_normalized, matched_points2_normalized);
  E = E_report.model;

  LORANSAC<FundamentalMatrixSevenPointEstimator,
           FundamentalMatrixEightPointEstimator>
      F_ransac(options.ransac_options);
  const auto F_report = F_ransac.Estimate(matched_points1, matched_points2);
  F = F_report.model;

  // Estimate planar or panoramic model.

  LORANSAC<HomographyMatrixEstimator, HomographyMatrixEstimator> H_ransac(
      options.ransac_options);
  const auto H_report = H_ransac.Estimate(matched_points1, matched_points2);
  H = H_report.model;

  if ((!E_report.success && !F_report.success && !H_report.success) ||
      (E_report.support.num_inliers < options.min_num_inliers &&
       F_report.support.num_inliers < options.min_num_inliers &&
       H_report.support.num_inliers < options.min_num_inliers)) {
    config = ConfigurationType::DEGENERATE;
    return;
  }

  // Determine inlier ratios of different models.

  const double E_F_inlier_ratio =
      static_cast<double>(E_report.support.num_inliers) /
      F_report.support.num_inliers;
  const double H_F_inlier_ratio =
      static_cast<double>(H_report.support.num_inliers) /
      F_report.support.num_inliers;
  const double H_E_inlier_ratio =
      static_cast<double>(H_report.support.num_inliers) /
      E_report.support.num_inliers;

  const std::vector<char>* best_inlier_mask = nullptr;
  size_t num_inliers = 0;

  if (E_report.success && E_F_inlier_ratio > options.min_E_F_inlier_ratio &&
      E_report.support.num_inliers >= options.min_num_inliers) {
    // Calibrated configuration.

    // Always use the model with maximum matches.
    if (E_report.support.num_inliers >= F_report.support.num_inliers) {
      num_inliers = E_report.support.num_inliers;
      best_inlier_mask = &E_report.inlier_mask;
    } else {
      num_inliers = F_report.support.num_inliers;
      best_inlier_mask = &F_report.inlier_mask;
    }

    if (H_E_inlier_ratio > options.max_H_inlier_ratio) {
      config = PLANAR_OR_PANORAMIC;
      if (H_report.support.num_inliers > num_inliers) {
        num_inliers = H_report.support.num_inliers;
        best_inlier_mask = &H_report.inlier_mask;
      }
    } else {
      config = ConfigurationType::CALIBRATED;
    }
  } else if (F_report.success &&
             F_report.support.num_inliers >= options.min_num_inliers) {
    // Uncalibrated configuration.

    num_inliers = F_report.support.num_inliers;
    best_inlier_mask = &F_report.inlier_mask;

    if (H_F_inlier_ratio > options.max_H_inlier_ratio) {
      config = ConfigurationType::PLANAR_OR_PANORAMIC;
      if (H_report.support.num_inliers > num_inliers) {
        num_inliers = H_report.support.num_inliers;
        best_inlier_mask = &H_report.inlier_mask;
      }
    } else {
      config = ConfigurationType::UNCALIBRATED;
    }
  } else if (H_report.success &&
             H_report.support.num_inliers >= options.min_num_inliers) {
    num_inliers = H_report.support.num_inliers;
    best_inlier_mask = &H_report.inlier_mask;
    config = ConfigurationType::PLANAR_OR_PANORAMIC;
  } else {
    config = ConfigurationType::DEGENERATE;
    return;
  }

  if (best_inlier_mask != nullptr) {
    inlier_matches =
        ExtractInlierMatches(matches, num_inliers, *best_inlier_mask);

    if (options.detect_watermark &&
        DetectWatermark(camera1, matched_points1, camera2, matched_points2,
                        num_inliers, *best_inlier_mask, options)) {
      config = ConfigurationType::WATERMARK;
    }
  }
}

// IS: GCRANSAC usage
//#pragma optimize("", off)
void TwoViewGeometry::EstimateUncalibratedGCRansac(
    const Camera& camera1, const std::vector<Eigen::Vector2d>& points1,
    const Camera& camera2, const std::vector<Eigen::Vector2d>& points2,
    const FeatureMatches& matches, const Options& options) {
  using namespace gcransac;

  options.Check();

  if (matches.size() < options.min_num_inliers) {
    config = ConfigurationType::DEGENERATE;
    return;
  }

  //GC RANSAC HERE
  // Algorithm Params
  const int fps = -1;  // The required FPS limit. If it is set to -1, the
                       // algorithm will not be interrupted before finishing.
  const double inlier_outlier_threshold_fundamental_matrix = 0.0005;  // The used adaptive inlier-outlier threshold in GC-RANSAC for
               // fundamental matrix estimation.
  const double inlier_outlier_threshold_homography = options.ransac_options.max_error;  // The used inlier-outlier threshold in GC-RANSAC for homography
             // estimation.
  const double spatial_coherence_weight = 0.14;  // The weight of the spatial coherence term in the graph-cut energy
             // minimization.
  const size_t cell_number_in_neighborhood_graph = 8;  // The number of cells along each axis in the neighborhood graph.

  cv::Mat cvPoints(static_cast<int>(matches.size()), 4, CV_64F);
  double* cvPointsPtr = reinterpret_cast<double*>(cvPoints.data);
  // IS: Transformation of the input
  // Extract corresponding points
  std::vector<Eigen::Vector2d> matched_points1(matches.size());
  std::vector<Eigen::Vector2d> matched_points2(matches.size());
  //std::vector<Eigen::Vector2d> matched_points1_normalized(matches.size());
  //std::vector<Eigen::Vector2d> matched_points2_normalized(matches.size());
  for (size_t i = 0; i < matches.size(); ++i) {
    const point2D_t idx1 = matches[i].point2D_idx1;
    const point2D_t idx2 = matches[i].point2D_idx2;
    *(cvPointsPtr++) = points1[idx1](0);
    *(cvPointsPtr++) = points1[idx1](1);
    *(cvPointsPtr++) = points2[idx2](0);
    *(cvPointsPtr++) = points2[idx2](1);
  }

  neighborhood::GridNeighborhoodGraph neighborhood(
      &cvPoints,
      camera1.Width() / static_cast<double>(cell_number_in_neighborhood_graph),
      camera1.Height() / static_cast<double>(cell_number_in_neighborhood_graph),
      camera2.Width() / static_cast<double>(cell_number_in_neighborhood_graph),
      camera2.Height() / static_cast<double>(cell_number_in_neighborhood_graph),
      cell_number_in_neighborhood_graph);

  // Checking if the neighborhood graph is initialized successfully.
  if (!neighborhood.isInitialized()) {
    fprintf(stderr,
            "The neighborhood graph is not initialized successfully.\n");
    return;
  }

  // Calculating the maximum image diagonal to be used for setting the threshold
  // adaptively for each image pair.
  const double max_image_diagonal =
      sqrt(pow(MAX(camera1.Width(), camera2.Width()), 2) +
           pow(MAX(camera1.Height(), camera2.Height()), 2));

  // Apply Graph-cut RANSAC for Fundamental Matrix
  utils::DefaultFundamentalMatrixEstimator estimator_F;
  //std::vector<int> inliers_F;
  FundamentalMatrix model_F;

  // Initialize the samplers
  // The main sampler is used inside the local optimization
  sampler::ProgressiveNapsacSampler main_sampler_F(
      &cvPoints,
      {16, 8, 4, 2},  // The layer of grids. The cells of the finest grid are of
                      // dimension (source_image_width / 16) *
                      // (source_image_height / 16)  * (destination_image_width
                      // / 16)  (destination_image_height / 16), etc.
      estimator_F.sampleSize(),                 // The size of a minimal sample
      static_cast<double>(camera1.Width()),   // The width of the source image
      static_cast<double>(camera1.Height()),  // The height of the source image
      static_cast<double>(camera2.Width()),  // The width of the destination image
      static_cast<double>(camera2.Height()));  // The height of the destination image

  // The local optimization sampler is used inside the local optimization
  sampler::UniformSampler local_optimization_sampler_F(&cvPoints);

  // Checking if the samplers are initialized successfully.
  if (!main_sampler_F.isInitialized() ||
      !local_optimization_sampler_F.isInitialized()) {
    fprintf(stderr, "One of the samplers is not initialized successfully.\n");
    return;
  }

  // 1. Fundamental Matrix
  GCRANSAC<utils::DefaultFundamentalMatrixEstimator,
           neighborhood::GridNeighborhoodGraph>
            gcransac_F;
  gcransac_F.setFPS(-1);  // Set the desired FPS (-1 means no limit)
  gcransac_F.settings.threshold = options.ransac_options.max_error;  // The inlier-outlier threshold
  gcransac_F.settings.spatial_coherence_weight = spatial_coherence_weight;  // The weight of the spatial coherence term
  gcransac_F.settings.confidence = options.ransac_options.confidence;  // The required confidence in the results
  gcransac_F.settings.max_local_optimization_number = 50;  // The maximum number of local optimizations
  gcransac_F.settings.max_iteration_number = options.ransac_options.max_num_trials;  // The maximum number of iterations
  gcransac_F.settings.min_iteration_number = options.ransac_options.min_num_trials;  // The minimum number of iterations
  gcransac_F.settings.neighborhood_sphere_radius = cell_number_in_neighborhood_graph;  // The radius of the neighborhood ball
  gcransac_F.settings.core_number = 1;  // The number of parallel processes

  // Start GC-RANSAC
  gcransac_F.run(cvPoints, 
               estimator_F, 
               &main_sampler_F, 
               &local_optimization_sampler_F,
               &neighborhood, model_F);

  // Get the statistics of the results
  const utils::RANSACStatistics& F_statistic = gcransac_F.getRansacStatistics();

  F = model_F.descriptor;

  //Homography GCRANSAC
  // Apply Graph-cut RANSAC for Essential Matrix
  utils::DefaultHomographyEstimator estimator_H;
  //std::vector<int> inliers;
  Homography model_H;

  // 1. Homography Matrix
  GCRANSAC<utils::DefaultHomographyEstimator, neighborhood::GridNeighborhoodGraph> gcransac_H;
  gcransac_H.setFPS(-1);  // Set the desired FPS (-1 means no limit)
  gcransac_H.settings.threshold = options.ransac_options.max_error;  // The inlier-outlier threshold
  gcransac_H.settings.spatial_coherence_weight = spatial_coherence_weight;  // The weight of the spatial coherence term
  gcransac_H.settings.confidence = options.ransac_options.confidence;  // The required confidence in the results
  gcransac_H.settings.max_local_optimization_number = 50;  // The maximum number of local optimizations
  gcransac_H.settings.max_iteration_number = options.ransac_options.max_num_trials;  // The maximum number of iterations
  gcransac_H.settings.min_iteration_number = options.ransac_options.min_num_trials;  // The minimum number of iterations
  gcransac_H.settings.neighborhood_sphere_radius = cell_number_in_neighborhood_graph;  // The radius of the neighborhood ball
  gcransac_H.settings.core_number = 1;    // The number of parallel processes

  // Initialize the samplers
  // The main sampler is used inside the local optimization
  sampler::ProgressiveNapsacSampler main_sampler_H(
      &cvPoints,
      {16, 8, 4, 2},  // The layer of grids. The cells of the finest grid are of
                      // dimension (source_image_width / 16) *
                      // (source_image_height / 16)  * (destination_image_width
                      // / 16)  (destination_image_height / 16), etc.
      estimator_H.sampleSize(),                 // The size of a minimal sample
      static_cast<double>(camera1.Width()),   // The width of the source image
      static_cast<double>(camera1.Height()),  // The height of the source image
      static_cast<double>(camera2.Width()),  // The width of the destination image
      static_cast<double>(camera2.Height()));  // The height of the destination image

  // The local optimization sampler is used inside the local optimization
  sampler::UniformSampler local_optimization_sampler_H(&cvPoints);

  // Checking if the samplers are initialized successfully.
  if (!main_sampler_H.isInitialized() ||
      !local_optimization_sampler_H.isInitialized()) {
    fprintf(stderr, "One of the samplers is not initialized successfully.\n");
    return;
  }

  // Start GC-RANSAC
  gcransac_H.run(cvPoints, 
               estimator_H, 
               &main_sampler_H, 
               &local_optimization_sampler_H,
               &neighborhood, 
               model_H);

  // Get the statistics of the results
  const utils::RANSACStatistics& H_statistic = gcransac_H.getRansacStatistics();

  H = model_H.descriptor;

  // TODO: fit statistic into report
  // if ((!F_report.success && !H_report.success) ||
  //    (F_report.support.num_inliers < options.min_num_inliers &&
  //     H_report.support.num_inliers < options.min_num_inliers)) {
  //  config = ConfigurationType::DEGENERATE;
  //  return;
  //}

  // Determine inlier ratios of different models.

  const double H_F_inlier_ratio =
      static_cast<double>(H_statistic.inliers.size()) /
      F_statistic.inliers.size();

  if (H_F_inlier_ratio > options.max_H_inlier_ratio) {
    config = ConfigurationType::PLANAR_OR_PANORAMIC;
  } else {
    config = ConfigurationType::UNCALIBRATED;
  }

  //IS: Create inlier_mask to get inliers in a way COLMAP has them
  std::vector<char> F_inlier_mask(matches.size(), false);
  for (auto idx : F_statistic.inliers) {
    F_inlier_mask[idx] = true;
  }

  inlier_matches = ExtractInlierMatches(matches, F_statistic.inliers.size(),
                                        F_inlier_mask);

  if (options.detect_watermark &&
      DetectWatermark(camera1, matched_points1, camera2, matched_points2,
                      F_statistic.inliers.size(), F_inlier_mask,
                      options)) {
    config = ConfigurationType::WATERMARK;
  }
}
//#pragma optimize("", on)


void TwoViewGeometry::EstimateUncalibrated(
    const Camera& camera1, const std::vector<Eigen::Vector2d>& points1,
    const Camera& camera2, const std::vector<Eigen::Vector2d>& points2,
    const FeatureMatches& matches, const Options& options) {
  options.Check();

  if (matches.size() < options.min_num_inliers) {
    config = ConfigurationType::DEGENERATE;
    return;
  }

  // Extract corresponding points.
  std::vector<Eigen::Vector2d> matched_points1(matches.size());
  std::vector<Eigen::Vector2d> matched_points2(matches.size());
  for (size_t i = 0; i < matches.size(); ++i) {
    matched_points1[i] = points1[matches[i].point2D_idx1];
    matched_points2[i] = points2[matches[i].point2D_idx2];
  }

  // Estimate epipolar model.

  LORANSAC<FundamentalMatrixSevenPointEstimator,
           FundamentalMatrixEightPointEstimator>
      F_ransac(options.ransac_options);
  const auto F_report = F_ransac.Estimate(matched_points1, matched_points2);
  F = F_report.model;

  // Estimate planar or panoramic model.

  LORANSAC<HomographyMatrixEstimator, HomographyMatrixEstimator> H_ransac(
      options.ransac_options);
  const auto H_report = H_ransac.Estimate(matched_points1, matched_points2);
  H = H_report.model;

  if ((!F_report.success && !H_report.success) ||
      (F_report.support.num_inliers < options.min_num_inliers &&
       H_report.support.num_inliers < options.min_num_inliers)) {
    config = ConfigurationType::DEGENERATE;
    return;
  }

  // Determine inlier ratios of different models.

  const double H_F_inlier_ratio =
      static_cast<double>(H_report.support.num_inliers) /
      F_report.support.num_inliers;

  if (H_F_inlier_ratio > options.max_H_inlier_ratio) {
    config = ConfigurationType::PLANAR_OR_PANORAMIC;
  } else {
    config = ConfigurationType::UNCALIBRATED;
  }

  inlier_matches = ExtractInlierMatches(matches, F_report.support.num_inliers,
                                        F_report.inlier_mask);

  if (options.detect_watermark &&
      DetectWatermark(camera1, matched_points1, camera2, matched_points2,
                      F_report.support.num_inliers, F_report.inlier_mask,
                      options)) {
    config = ConfigurationType::WATERMARK;
  }
}

bool TwoViewGeometry::DetectWatermark(
    const Camera& camera1, const std::vector<Eigen::Vector2d>& points1,
    const Camera& camera2, const std::vector<Eigen::Vector2d>& points2,
    const size_t num_inliers, const std::vector<char>& inlier_mask,
    const Options& options) {
  options.Check();

  // Check if inlier points in border region and extract inlier matches.

  const double diagonal1 = std::sqrt(camera1.Width() * camera1.Width() +
                                     camera1.Height() * camera1.Height());
  const double diagonal2 = std::sqrt(camera2.Width() * camera2.Width() +
                                     camera2.Height() * camera2.Height());
  const double minx1 = options.watermark_border_size * diagonal1;
  const double miny1 = minx1;
  const double maxx1 = camera1.Width() - minx1;
  const double maxy1 = camera1.Height() - miny1;
  const double minx2 = options.watermark_border_size * diagonal2;
  const double miny2 = minx2;
  const double maxx2 = camera2.Width() - minx2;
  const double maxy2 = camera2.Height() - miny2;

  std::vector<Eigen::Vector2d> inlier_points1(num_inliers);
  std::vector<Eigen::Vector2d> inlier_points2(num_inliers);

  size_t num_matches_in_border = 0;

  size_t j = 0;
  for (size_t i = 0; i < inlier_mask.size(); ++i) {
    if (inlier_mask[i]) {
      const auto& point1 = points1[i];
      const auto& point2 = points2[i];

      inlier_points1[j] = point1;
      inlier_points2[j] = point2;
      j += 1;

      if (!IsImagePointInBoundingBox(point1, minx1, maxx1, miny1, maxy1) &&
          !IsImagePointInBoundingBox(point2, minx2, maxx2, miny2, maxy2)) {
        num_matches_in_border += 1;
      }
    }
  }

  const double matches_in_border_ratio =
      static_cast<double>(num_matches_in_border) / num_inliers;

  if (matches_in_border_ratio < options.watermark_min_inlier_ratio) {
    return false;
  }

  // Check if matches follow a translational model.

  RANSACOptions ransac_options = options.ransac_options;
  ransac_options.min_inlier_ratio = options.watermark_min_inlier_ratio;

  LORANSAC<TranslationTransformEstimator<2>, TranslationTransformEstimator<2>>
      ransac(ransac_options);
  const auto report = ransac.Estimate(inlier_points1, inlier_points2);

  const double inlier_ratio =
      static_cast<double>(report.support.num_inliers) / num_inliers;

  return inlier_ratio >= options.watermark_min_inlier_ratio;
}

}  // namespace colmap
