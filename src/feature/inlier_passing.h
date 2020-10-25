//Ilia Shipachev, CVUT-FEL, CMP lab, 2020

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

  void save_inliers(image_t img1, image_t img2, FeatureMatches &inliers);
  void reorder_by_passed_inliers(image_t img_j,
                                 image_t img_k,
                                 FeatureMatches& matches);
  //get list of amount of inliers passed to img_j from all other img_i for all i's
  std::unordered_map<image_t, size_t> const calc_inliers_passed(image_t img_j, image_t img_k);
  std::unordered_map<image_t, std::vector<size_t>> const get_inliers_passed(image_t img_j, image_t img_k);

 private:
   //j's inliers of model i-to-j related to j's detections
  std::unordered_map<std::pair<image_t, image_t>, std::vector<point2D_t>> pair_inliers; 
  std::vector<std::vector<image_t>> connected_by; //list of all models builded to i-s image

  typedef std::tuple<image_t, image_t, image_t> triplet_t;

  //size_t hash_fn (const triplet_t& t) {
  //  size_t res = 17;
  //  res = res * 31 + std::hash<image_t>()(std::get<0>(t));
  //  res = res * 31 + std::hash<image_t>()(std::get<1>(t));
  //  res = res * 31 + std::hash<image_t>()(std::get<2>(t));
  //  return res;  
  //};

  struct hash_fn : public std::unary_function<triplet_t, size_t> {
    size_t operator()(const triplet_t& t) const {
      size_t res = 17;
      res = res * 31 + std::hash<image_t>()(std::get<0>(t));
      res = res * 31 + std::hash<image_t>()(std::get<1>(t));
      res = res * 31 + std::hash<image_t>()(std::get<2>(t));
      return res;
    };
  };


  std::unordered_map<triplet_t, std::vector<size_t>, hash_fn> pci;
};

}  // namespace colmap
