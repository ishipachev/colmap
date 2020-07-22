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

// IS: MAGSAC additions
#include "magsac.h"
#include "magsac/graph-cut-ransac/src/pygcransac/include/utils.h"
#include "magsac/include/estimators.h"
#include <opencv2/core/core.hpp>

#pragma message("warning: Hi")

// IS: GCRansac additions
//#include "GCRANSAC.h"
//#include "essential_estimator.h"
//#include "fundamental_estimator.h"
//#include "gcransac_utils.h"
//#include "grid_neighborhood_graph.h"
//#include "homography_estimator.h"
//#include "model.h"
//#include "progressive_napsac_sampler.h"
//#include "uniform_sampler.h"
//#include <opencv2/core/core.hpp>

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

bool setAndRunMAGSAC_F(const cv::Mat& cvPoints,
                       const TwoViewGeometry::Options& options,
                       Eigen::Matrix3d& F, std::vector<char>& inlier_mask,
                       size_t& inlier_cnt) {
  // Estimate epipolar model.
  // MAGSAC way
  gcransac::sampler::UniformSampler main_sampler_F(&cvPoints);
  gcransac::FundamentalMatrix model_F;
  ModelScore score;
  int iteration_number = 0;  // Number of iterations required

  MAGSAC<cv::Mat, magsac::utils::DefaultFundamentalMatrixEstimator> magsac_F;
  // max_sigma will be twice max error
  magsac::utils::DefaultFundamentalMatrixEstimator estimator_F(
      options.ransac_options.max_error * 4);
  magsac_F.setMaximumThreshold(options.ransac_options.max_error * 4);
  magsac_F.setReferenceThreshold(options.ransac_options.max_error);
  magsac_F.setIterationLimit(options.ransac_options.max_num_trials);
  magsac_F.setMinimumIterationNumber(options.ransac_options.min_num_trials);
  magsac_F.setFPS(0);

  bool success =
      magsac_F.run(cvPoints,  // The data points
                   options.ransac_options.confidence,   // The required confidence in the results
                   estimator_F,       // The used estimator
                   main_sampler_F,    // The sampler used for selecting minimal
                                      // samples in each iteration
                   model_F,           // The estimated model
                   iteration_number,  // The number of iterations
                   score);            // The score of the estimated model

  F = model_F.descriptor;

  std::vector<bool> inlier_mask_bool(cvPoints.rows);

  magsac_F.getModelInliersMask(cvPoints, model_F, estimator_F,
                               options.ransac_options.max_error,
                               inlier_mask_bool);
  inlier_cnt = static_cast<uint>(
							   std::count(inlier_mask_bool.begin(), inlier_mask_bool.end(), true));
  if (inlier_mask.size() != cvPoints.rows) {
    inlier_mask.resize(cvPoints.rows);
  }

  // converting to chars
  for (uint i = 0; i < cvPoints.rows; ++i) {
    inlier_mask[i] = static_cast<char>(inlier_mask_bool[i]);
  }

  return success;
}

bool setAndRunMAGSAC_H(const cv::Mat& cvPoints,
                       const TwoViewGeometry::Options& options,
                       Eigen::Matrix3d& H, std::vector<char>& inlier_mask,
                       size_t& inlier_cnt) {
  // Estimate epipolar model.
  // MAGSAC way
  gcransac::sampler::UniformSampler main_sampler_F(&cvPoints);
  gcransac::Homography model_H;
  ModelScore score;
  magsac::utils::DefaultHomographyEstimator estimator_H;
  MAGSAC<cv::Mat, magsac::utils::DefaultHomographyEstimator> magsac_H;
  // max_sigma will be twice max error
  magsac_H.setMaximumThreshold(options.ransac_options.max_error * 4);
  magsac_H.setReferenceThreshold(options.ransac_options.max_error);
  magsac_H.setIterationLimit(options.ransac_options.max_num_trials);
  magsac_H.setMinimumIterationNumber(options.ransac_options.min_num_trials);
  magsac_H.setFPS(0);

  int iteration_number = 0;  // Number of iterations required

  bool success =
      magsac_H.run(cvPoints,  // The data points
                   options.ransac_options.confidence,   // The required confidence in the results
                   estimator_H,       // The used estimator
                   main_sampler_F,    // The sampler used for selecting minimal
                                      // samples in each iteration
                   model_H,           // The estimated model
                   iteration_number,  // The number of iterations
                   score);            // The score of the estimated model

  H = model_H.descriptor;

  std::vector<bool> inlier_mask_bool(cvPoints.rows);

  magsac_H.getModelInliersMask(cvPoints, model_H, estimator_H,
                               options.ransac_options.max_error,
                               inlier_mask_bool);
  inlier_cnt = static_cast<uint>(
      std::count(inlier_mask_bool.begin(), inlier_mask_bool.end(), true));

  if (inlier_mask.size() != cvPoints.rows) {
    inlier_mask.resize(cvPoints.rows);
  }

  // converting to chars
  for (uint i = 0; i < cvPoints.rows; ++i) {
    inlier_mask[i] = static_cast<char>(inlier_mask_bool[i]);
  }

  return success;
}

}  // namespace

void TwoViewGeometry::EstimateMAGSAC(const Camera& camera1,
                               const std::vector<Eigen::Vector2d>& points1,
                               const Camera& camera2,
                               const std::vector<Eigen::Vector2d>& points2,
                               const FeatureMatches& matches,
                               const Options& options) {
  if (camera1.HasPriorFocalLength() && camera2.HasPriorFocalLength()) {
    EstimateCalibratedMAGSAC(camera1, points1, camera2, points2, matches, options);
  } else {
    EstimateUncalibratedMAGSAC(camera1, points1, camera2, points2, matches, options);
  }
}

void TwoViewGeometry::EstimateMultipleMAGSAC(
    const Camera& camera1, const std::vector<Eigen::Vector2d>& points1,
    const Camera& camera2, const std::vector<Eigen::Vector2d>& points2,
    const FeatureMatches& matches, const Options& options) {
  FeatureMatches remaining_matches = matches;
  std::vector<TwoViewGeometry> two_view_geometries;
  while (true) {
    TwoViewGeometry two_view_geometry;
    two_view_geometry.EstimateMAGSAC(camera1, points1, camera2, points2,
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

void TwoViewGeometry::EstimateCalibratedMAGSAC(
    const Camera& camera1, const std::vector<Eigen::Vector2d>& points1,
    const Camera& camera2, const std::vector<Eigen::Vector2d>& points2,
    const FeatureMatches& matches, const Options& options) {
  options.Check();

  printf("\t EstimateCalibratedMAGSAC: function is not done yet!\n");

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
        DetectWatermarkMAGSAC(camera1, matched_points1, camera2,
                              matched_points2, num_inliers, *best_inlier_mask,
                              options)) {
      config = ConfigurationType::WATERMARK;
    }
  }
}

void TwoViewGeometry::EstimateUncalibratedMAGSAC(
    const Camera& camera1, const std::vector<Eigen::Vector2d>& points1,
    const Camera& camera2, const std::vector<Eigen::Vector2d>& points2,
    const FeatureMatches& matches, const Options& options) {
  options.Check();

  if (matches.size() < options.min_num_inliers) {
    config = ConfigurationType::DEGENERATE;
    return;
  }

  //printf("\tUncalibrated MAGSAC is running...\n");

  // Extract corresponding points.
  // And transforming the input
  cv::Mat cvPoints(static_cast<int>(matches.size()), 4, CV_64F);
  double* cvPointsPtr = reinterpret_cast<double*>(cvPoints.data);

  std::vector<Eigen::Vector2d> matched_points1(matches.size());
  std::vector<Eigen::Vector2d> matched_points2(matches.size());
  for (size_t i = 0; i < matches.size(); ++i) {
    matched_points1[i] = points1[matches[i].point2D_idx1];
    matched_points2[i] = points2[matches[i].point2D_idx2];

    const point2D_t idx1 = matches[i].point2D_idx1;
    const point2D_t idx2 = matches[i].point2D_idx2;
    *(cvPointsPtr++) = points1[idx1](0);
    *(cvPointsPtr++) = points1[idx1](1);
    *(cvPointsPtr++) = points2[idx2](0);
    *(cvPointsPtr++) = points2[idx2](1);
  }

  // Estimate epipolar model
  // Default way
  LORANSAC<FundamentalMatrixSevenPointEstimator,
           FundamentalMatrixEightPointEstimator>
      F_ransac(options.ransac_options);
  const auto F_report = F_ransac.Estimate(matched_points1, matched_points2);
  F = F_report.model;

  // Estimate epipolar model.
  // MAGSAC way
  // Using default ransac::report structure to have a report
  LORANSAC<FundamentalMatrixSevenPointEstimator,
           FundamentalMatrixEightPointEstimator>::Report F_magsac_report;
  F_magsac_report.success = setAndRunMAGSAC_F(cvPoints, 
                                              options, 
                                              F, 
                                              F_magsac_report.inlier_mask, 
                                              F_magsac_report.support.num_inliers);

  // Estimate planar or panoramic model.
  // Default way
  LORANSAC<HomographyMatrixEstimator, HomographyMatrixEstimator> H_ransac(
      options.ransac_options);
  const auto H_report = H_ransac.Estimate(matched_points1, matched_points2);
  H = H_report.model;

  // Estimate planar or panoramic model.
  // MAGSAC way
  // Using default ransac::report structure to have a report
  LORANSAC<HomographyMatrixEstimator, 
           HomographyMatrixEstimator>::Report H_magsac_report;
  H_magsac_report.success = setAndRunMAGSAC_H(cvPoints, 
                                              options, 
                                              H, 
                                              H_magsac_report.inlier_mask, 
                                              H_magsac_report.support.num_inliers);

  if ((!F_magsac_report.success && !H_magsac_report.success) ||
      (F_magsac_report.support.num_inliers < options.min_num_inliers &&
       H_magsac_report.support.num_inliers < options.min_num_inliers)) {
    config = ConfigurationType::DEGENERATE;
    return;
  }

  printf("F, default   : %5d\n", F_report.support.num_inliers);
  printf("F, magsac    : %5d\n", F_magsac_report.support.num_inliers);
  printf("H, default: %5d\n", H_report.support.num_inliers);
  printf("H, magsac : %5d\n", H_magsac_report.support.num_inliers);


  // Determine inlier ratios of different models.

  const double H_F_inlier_ratio =
      static_cast<double>(H_magsac_report.support.num_inliers) /
      F_magsac_report.support.num_inliers;

  if (H_F_inlier_ratio > options.max_H_inlier_ratio) {
    config = ConfigurationType::PLANAR_OR_PANORAMIC;
  } else {
    config = ConfigurationType::UNCALIBRATED;
  }

  inlier_matches = ExtractInlierMatches(
      matches, 
      F_magsac_report.support.num_inliers, 
      F_magsac_report.inlier_mask);

  if (options.detect_watermark &&
      DetectWatermarkMAGSAC(camera1, matched_points1, camera2, matched_points2,
                            F_magsac_report.support.num_inliers,
                            F_magsac_report.inlier_mask,
                            options)) {
    config = ConfigurationType::WATERMARK;
  }
}

bool TwoViewGeometry::DetectWatermarkMAGSAC(
    const Camera& camera1, const std::vector<Eigen::Vector2d>& points1,
    const Camera& camera2, const std::vector<Eigen::Vector2d>& points2,
    const size_t num_inliers, const std::vector<char>& inlier_mask,
    const Options& options) {
  options.Check();

  printf("Detecting watermark (MAGSAC is not supported here\n");

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
