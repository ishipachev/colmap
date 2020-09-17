#pragma once
#include "util/types.h"
#include "feature/types.h"
#include <unordered_map>
#include <vector>

#define INL_PASSING_ON false

namespace colmap {

class InlierPassing {
 public:
  InlierPassing();

  void write_inliers(image_t img1, image_t img2, FeatureMatches &inliers);
  void reorder_matches_by_passed_inliers(image_t img_j,
                                         image_t img_k,
                                         FeatureMatches& matches);

 private:
  //std::unordered_map<std::pair<image_t, image_t>, FeatureMatches> pair_inliers; //inliers of model i-to-j
  std::unordered_map<std::pair<image_t, image_t>, std::vector<point2D_t>> pair_inliers; //j's inliers of model i-to-j

  std::vector<std::vector<image_t>> connected_by; //list of all models builded to i-s image
};

}  // namespace colmap