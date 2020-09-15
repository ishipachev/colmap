#include "inlier_passing.h"
#include "feature/types.h"
#include "util/logging.h"

#pragma optimize( "", off )

namespace colmap {

//calculate intersection indicies related to matches array
//both 1d arrays inliers and matches are expected to be sorted
void get_intersection_ids(std::vector<point2D_t> &inliers, 
                          std::vector<point2D_t> &matches, 
                          std::vector<size_t> &pci) {
  size_t i_idx = 0;
  size_t m_idx = 0;

  point2D_t i;
  point2D_t m;

  pci.resize(inliers.size());
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
}

InlierPassing::InlierPassing(){
}

void InlierPassing::write_inliers(image_t img_i, image_t img_j, FeatureMatches& inliers) {

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
void InlierPassing::reorder_matches_by_passed_inliers(image_t img_j, 
                                                      image_t img_k, 
                                                      FeatureMatches& matches_jk) {
  printf("\nReordering matches for images %d and %d\n", img_j, img_k);
  //get rid of second index which we don't need
  //suppose to be sorted by point2D_idx1 field
  std::vector<point2D_t> matches_jk_j(matches_jk.size());
  for (int s = 0; s < matches_jk_j.size(); ++s) {
    matches_jk_j[s] = matches_jk[s].point2D_idx1;   //from pair copying only indicies of j's image 
  }
  CHECK(std::is_sorted(matches_jk_j.begin(), matches_jk_j.end()));

  if (connected_by.size() > img_j) {
    size_t reord_idx = 0;
    for (image_t img_i : connected_by[img_j]) {      
      std::vector<point2D_t> &inliers_ij_j = pair_inliers[{img_i, img_j}];
      std::vector<size_t> pci_ijk;
      get_intersection_ids(inliers_ij_j, matches_jk_j, pci_ijk);
      for (auto s: pci_ijk) {
        std::swap(matches_jk[reord_idx++], matches_jk[s]);
      }
      printf("%d matches were reordered\n using model from %d to %d\n", inl_idx, img_i, img_j);
    }
  }
  else {
    printf("No existing models to image %d", img_j);
  }
}

//bool InlierPassing::is_geom_found(image_t img1, image_t img2) { 
//  return static_cast<bool>(pair_inliers.count({img1, img2})); 
//};




} // namespace colmap