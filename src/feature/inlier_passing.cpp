//Ilia Shipachev, CVUT-FEL, CMP lab, 2020

#include "inlier_passing.h"
#include "feature/types.h"
#include "util/logging.h"
#include <vector>

//#pragma optimize( "", off )

namespace colmap {

//calculate intersection indicies related to matches array
//both 1d arrays inliers and matches are expected to be sorted
//const std::vector<size_t> &get_intersection_ids(std::vector<point2D_t> &inliers, 
std::vector<size_t> get_intersection_ids(std::vector<point2D_t> &inliers, 
                                         std::vector<point2D_t> &matches) {
  size_t i_idx = 0;
  size_t m_idx = 0;

  point2D_t i;
  point2D_t m;

  std::vector<size_t> pci(inliers.size());
  size_t pci_idx = 0;

  while ((i_idx < inliers.size()) && (m_idx < inliers.size())) {
    i = inliers[i_idx];
    m = matches[m_idx];

    if (i == m) {
      pci[pci_idx++] = m_idx;
      i_idx++;
      m_idx++;
    } else if (i > m) {
      m_idx++;
    } else {
      i_idx++;
    }
  }

  pci.resize(pci_idx);
  return pci;
}

InlierPassing::InlierPassing(){
}

void InlierPassing::save_inliers(image_t img_i, image_t img_j, FeatureMatches& inliers) {

  //pair_inliers[{img_i, img_j}].resize(inliers.size());

  std::vector<point2D_t> &pair_inliers_ij = pair_inliers[{img_i, img_j}];
  pair_inliers_ij.resize(inliers.size());

  for (int s = 0; s < pair_inliers_ij.size(); ++s) {
    pair_inliers_ij[s] = inliers[s].point2D_idx2;   //from pair copying only indicies of j's image 
  }
  std::sort(pair_inliers_ij.begin(), pair_inliers_ij.end());

  if (connected_by.size() < img_j + 1) { connected_by.resize(img_j + 1); }
  connected_by[img_j].push_back(img_i);
}

//Reordering matches_jk by putting first inliers from aldready calculated models from i to j
//for all such i
void InlierPassing::reorder_by_passed_inliers(image_t img_j, 
                                              image_t img_k, 
                                              FeatureMatches& matches_jk) {

//TODO: Rewrite into: check the size, is smaller -- resize
//if bigger, at fist check the connected_by[img_j].size()
//otherwise no point to sort matches and do other unnececarliy job

  if (connected_by.size() > img_j) {
    printf("Reordering matches for images %d and %d\n", img_j, img_k);
    //get rid of the second index which we don't need
    //suppose to be sorted by point2D_idx1 field
    std::vector<point2D_t> matches_jk_j(matches_jk.size());
    for (int s = 0; s < matches_jk_j.size(); ++s) {
      matches_jk_j[s] = matches_jk[s].point2D_idx1;   //from pair copying only indicies of j's image 
    }
    CHECK(std::is_sorted(matches_jk_j.begin(), matches_jk_j.end()));

    size_t reord_idx = 0;
    for (image_t img_i : connected_by[img_j]) {      
      std::vector<point2D_t> &inliers_ij_j = pair_inliers[{img_i, img_j}];
      std::vector<size_t> pci_ijk;
      //std::vector<size_t> pci_ijk = get_intersection_ids(inliers_ij_j, matches_jk_j, pci_ijk);
      pci_ijk = get_intersection_ids(inliers_ij_j, matches_jk_j);
      pci[{img_i, img_j, img_k}] = pci_ijk;
      for (auto s: pci_ijk) {
        std::swap(matches_jk[reord_idx++], matches_jk[s]);
      }
      printf("%d matches were reordered of %d total matches by using model from %d to %d\n", reord_idx, matches_jk.size(), img_i, img_j);
      if (reord_idx >= matches_jk.size()) {
        break; //this part should be rewritten that prioritize points which we often have as inliers for other models
               //maybe just implement as a simple counter and order them in this way
      }
      break; //right now we will do only one iteratin for the testing purpose
    }
  } else {
    connected_by.resize(img_j + 1);
    printf("No existing models to image %d\n", img_j);
  }
}

std::unordered_map<image_t, size_t> 
const InlierPassing::calc_inliers_passed(image_t img_j, image_t img_k) {
  std::unordered_map<image_t, size_t> res;
  if (connected_by.size() >= img_j + 1) {
    for (auto img_i : connected_by[img_j]) {
      CHECK(img_i < img_j);   //not implemented for img_i > img_j
      CHECK(img_j < img_k);   //not implemented for img_j > img_k
      res[img_i] = pci[{img_i, img_j, img_k}].size();
    }
  }
  return res;
}

std::unordered_map<image_t, std::vector<size_t>> 
const InlierPassing::get_inliers_passed(image_t img_j, image_t img_k) {
  std::unordered_map<image_t, std::vector<size_t>> res;
  if (connected_by.size() >= img_j + 1) {
    for (auto img_i : connected_by[img_j]) {
      CHECK(img_i < img_j);   //not implemented for img_i > img_j
      CHECK(img_j < img_k);   //not implemented for img_j > img_k
      res[img_i] = pci[{img_i, img_j, img_k}];
    }
  }
  return res;
}

//bool InlierPassing::is_geom_found(image_t img1, image_t img2) { 
//  return static_cast<bool>(pair_inliers.count({img1, img2})); 
//};




} // namespace colmap