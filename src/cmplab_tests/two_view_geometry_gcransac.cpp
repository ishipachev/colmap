// Copyright (C) 2019 Czech Technical University
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of Czech Technical University nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
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
// Please contact the author of this library if you have any questions.
// Author: Ilia Shipachev (shipaili@fel.cvut.cz)
// Date: 2019.11.01, CMP Lab

#include "base/camera.h"
#include "estimators/two_view_geometry.h"
#include "feature/types.h"
#include "optim/ransac.h"

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

//--------HELPER FUNCTIONS-----------//
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
}  // namespace

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
  const double spatial_coherence_weight = 0.14;  // The weight of the spatial coherence term in the graph-cut energy
             // minimization.
  const size_t cell_number_in_neighborhood_graph = 8;  // The number of cells along each axis in the neighborhood graph.

  cv::Mat cvPoints(static_cast<int>(matches.size()), 4, CV_64F);
  double* cvPointsPtr = reinterpret_cast<double*>(cvPoints.data);
  // IS: Transformation of the input
  // Extract corresponding points
  std::vector<Eigen::Vector2d> matched_points1(matches.size());
  std::vector<Eigen::Vector2d> matched_points2(matches.size());
  for (size_t i = 0; i < matches.size(); ++i) {
    const point2D_t idx1 = matches[i].point2D_idx1;
    const point2D_t idx2 = matches[i].point2D_idx2;
    *(cvPointsPtr++) = points1[idx1](0);
    *(cvPointsPtr++) = points1[idx1](1);
    *(cvPointsPtr++) = points2[idx2](0);
    *(cvPointsPtr++) = points2[idx2](1);
  }

  neighborhood::GridNeighborhoodGraph neighborhood_F(
      &cvPoints,
      camera1.Width() / static_cast<double>(cell_number_in_neighborhood_graph),
      camera1.Height() / static_cast<double>(cell_number_in_neighborhood_graph),
      camera2.Width() / static_cast<double>(cell_number_in_neighborhood_graph),
      camera2.Height() / static_cast<double>(cell_number_in_neighborhood_graph),
      cell_number_in_neighborhood_graph);

  // Checking if the neighborhood graph is initialized successfully.
  if (!neighborhood_F.isInitialized()) {
    fprintf(stderr,
            "The neighborhood graph is not initialized successfully for F.\n");
    return;
  }

  // Apply Graph-cut RANSAC for Fundamental Matrix
  utils::DefaultFundamentalMatrixEstimator estimator_F;
  //std::vector<int> inliers_F;
  FundamentalMatrix model_F;

  // Initialize the samplers
  // The main sampler is used inside the local optimization
  //sampler::ProgressiveNapsacSampler main_sampler_F(
  //    &cvPoints,
  //    {16, 8, 4, 2},  // The layer of grids. The cells of the finest grid are of
  //                    // dimension (source_image_width / 16) *
  //                    // (source_image_height / 16)  * (destination_image_width
  //                    // / 16)  (destination_image_height / 16), etc.
  //    estimator_F.sampleSize(),                 // The size of a minimal sample
  //    static_cast<double>(camera1.Width()),   // The width of the source image
  //    static_cast<double>(camera1.Height()),  // The height of the source image
  //    static_cast<double>(camera2.Width()),  // The width of the destination image
  //    static_cast<double>(camera2.Height()));  // The height of the destination image

  sampler::UniformSampler main_sampler_F(&cvPoints);
  // The local optimization sampler is used inside the local optimization
  sampler::UniformSampler local_optimization_sampler_F(&cvPoints);

  // Checking if the samplers are initialized successfully.
  if (!main_sampler_F.isInitialized() ||
      !local_optimization_sampler_F.isInitialized()) {
    fprintf(stderr, "One of the samplers is not initialized successfully for F.\n");
    return;
  }

  // 1. Fundamental Matrix
  GCRANSAC<utils::DefaultFundamentalMatrixEstimator,
           neighborhood::GridNeighborhoodGraph>
            gcransac_F;
  gcransac_F.setFPS(-1);  // Set the desired FPS (-1 means no limit)
  gcransac_F.settings.threshold = sqrt(options.ransac_options.max_error);  // The inlier-outlier threshold
  gcransac_F.settings.spatial_coherence_weight = spatial_coherence_weight;  // The weight of the spatial coherence term
  gcransac_F.settings.confidence = options.ransac_options.confidence;  // The required confidence in the results
  gcransac_F.settings.max_local_optimization_number = 1000;  // The maximum number of local optimizations
  gcransac_F.settings.max_iteration_number = options.ransac_options.max_num_trials;  // The maximum number of iterations
  gcransac_F.settings.min_iteration_number = options.ransac_options.min_num_trials;  // The minimum number of iterations
  gcransac_F.settings.neighborhood_sphere_radius = cell_number_in_neighborhood_graph;  // The radius of the neighborhood ball
  gcransac_F.settings.core_number = 1;  // The number of parallel processes

  // Start GC-RANSAC
  gcransac_F.run(cvPoints, 
               estimator_F, 
               &main_sampler_F, 
               &local_optimization_sampler_F,
               &neighborhood_F, model_F);

  // Get the statistics of the results
  const utils::RANSACStatistics& F_statistic = gcransac_F.getRansacStatistics();

  F = model_F.descriptor;

  //-------------------------------
  //Homography GCRANSAC
  //-------------------------------

  neighborhood::GridNeighborhoodGraph neighborhood_H(
    &cvPoints,
    camera1.Width() / static_cast<double>(cell_number_in_neighborhood_graph),
    camera1.Height() / static_cast<double>(cell_number_in_neighborhood_graph),
    camera2.Width() / static_cast<double>(cell_number_in_neighborhood_graph),
    camera2.Height() / static_cast<double>(cell_number_in_neighborhood_graph),
    cell_number_in_neighborhood_graph);

  // Checking if the neighborhood graph is initialized successfully.
  if (!neighborhood_H.isInitialized()) {
    fprintf(stderr,
            "The neighborhood graph is not initialized successfully for H.\n");
    return;
  }

  // Apply Graph-cut RANSAC for Essential Matrix
  utils::DefaultHomographyEstimator estimator_H;
  //std::vector<int> inliers;
  Homography model_H;

  // 1. Homography Matrix
  GCRANSAC<utils::DefaultHomographyEstimator, neighborhood::GridNeighborhoodGraph> gcransac_H;
  gcransac_H.setFPS(-1);  // Set the desired FPS (-1 means no limit)
  gcransac_H.settings.threshold = sqrt(options.ransac_options.max_error);  // The inlier-outlier threshold
  gcransac_H.settings.spatial_coherence_weight = spatial_coherence_weight;  // The weight of the spatial coherence term
  gcransac_H.settings.confidence = options.ransac_options.confidence;  // The required confidence in the results
  gcransac_H.settings.max_local_optimization_number = 1000;  // The maximum number of local optimizations
  gcransac_H.settings.max_iteration_number = options.ransac_options.max_num_trials;  // The maximum number of iterations
  gcransac_H.settings.min_iteration_number = options.ransac_options.min_num_trials;  // The minimum number of iterations
  gcransac_H.settings.neighborhood_sphere_radius = cell_number_in_neighborhood_graph;  // The radius of the neighborhood ball
  gcransac_H.settings.core_number = 1;    // The number of parallel processes

  // Initialize the samplers
  // The main sampler is used inside the local optimization
  //sampler::ProgressiveNapsacSampler main_sampler_H(
  //    &cvPoints,
  //    {16, 8, 4, 2},  // The layer of grids. The cells of the finest grid are of
  //                    // dimension (source_image_width / 16) *
  //                    // (source_image_height / 16)  * (destination_image_width
  //                    // / 16)  (destination_image_height / 16), etc.
  //    estimator_H.sampleSize(),                 // The size of a minimal sample
  //    static_cast<double>(camera1.Width()),   // The width of the source image
  //    static_cast<double>(camera1.Height()),  // The height of the source image
  //    static_cast<double>(camera2.Width()),  // The width of the destination image
  //    static_cast<double>(camera2.Height()), // The height of the destination image
  //    0.5);  // The length (i.e., 0.5 * <point number> iterations) of fully blending to global sampling 

  sampler::UniformSampler main_sampler_H(&cvPoints);
  // The local optimization sampler is used inside the local optimization
  sampler::UniformSampler local_optimization_sampler_H(&cvPoints);

  // Checking if the samplers are initialized successfully.
  if (!main_sampler_H.isInitialized() ||
      !local_optimization_sampler_H.isInitialized()) {
    fprintf(stderr, "One of the samplers is not initialized successfully for H.\n");
    return;
  }

  // Start GC-RANSAC
  gcransac_H.run(cvPoints, 
               estimator_H, 
               &main_sampler_H, 
               &local_optimization_sampler_H,
               &neighborhood_H, 
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



}  // namespace colmap