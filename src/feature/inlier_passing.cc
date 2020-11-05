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

  while ((i_idx < inliers.size()) && (m_idx < matches.size())) {
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
  std::vector<point2D_t> &pair_inliers_ji = pair_inliers[{img_j, img_i}];

  pair_inliers_ij.resize(inliers.size());
  pair_inliers_ji.resize(inliers.size());

  for (size_t s = 0; s < pair_inliers_ij.size(); ++s) {
    pair_inliers_ij[s] = inliers[s].point2D_idx2;   //from pair copying only indicies of j's image 
    pair_inliers_ji[s] = inliers[s].point2D_idx1;   //from pair copying only indiceis of i's image
  }
  std::sort(pair_inliers_ij.begin(), pair_inliers_ij.end());
  std::sort(pair_inliers_ji.begin(), pair_inliers_ji.end());

  if (connected_by.size() < img_j + 1) { connected_by.resize(img_j + 1); }
  if (connected_by.size() < img_i + 1) { connected_by.resize(img_i + 1); }

  connected_by[img_j].push_back(img_i);
  connected_by[img_i].push_back(img_j);
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
    //supposed to be sorted by point2D_idx1 field
    std::vector<point2D_t> matches_jk_j(matches_jk.size());
    for (size_t s = 0; s < matches_jk_j.size(); ++s) {
      matches_jk_j[s] = matches_jk[s].point2D_idx1;   //from pair copying only indicies of j's image 
    }
    //so just by all experiments its always sorted as it should be
    CHECK(std::is_sorted(matches_jk_j.begin(), matches_jk_j.end()));

    size_t best_i = 0;
    //size_t best_pci_size = 0;
    std::vector<size_t> best_pci_ijk(0);
    for (image_t img_i : connected_by[img_j]) {      
      std::vector<size_t> pci_ijk;

      std::vector<point2D_t> &inliers_ij_j = pair_inliers[{img_i, img_j}];
      if (inliers_ij_j.size() <= best_pci_ijk.size()) {
        //that means that with intersection we won't get more than in best
	//as it is checked further, so just skip it to save time
	continue; 
      }
      CHECK(std::is_sorted(inliers_ij_j.begin(), inliers_ij_j.end()));
      pci_ijk = get_intersection_ids(inliers_ij_j, matches_jk_j);

      CHECK(std::is_sorted(pci_ijk.begin(), pci_ijk.end()));
      //substituted with just a sort check, because it should be sorted already
      //std::sort(pci_ijk.begin(), pci_ijk.end()); 

      //pci[{img_i, img_j, img_k}] = pci_ijk;
      //if (best_pci_size < pci_ijk.size()){
      //  best_pci_size = pci_ijk.size();
      //  best_i = img_i;
      if (best_pci_ijk.size() < pci_ijk.size()){
	best_pci_ijk = std::move(pci_ijk);
	best_i = img_i;
      }
    }      
    pci[{best_i, img_j, img_k}] = best_pci_ijk;
    size_t reord_idx = 0;
    for (auto s: pci[{best_i, img_j, img_k}]) {
      if (reord_idx >= matches_jk.size() || s >= matches_jk.size()) { //Not needed but just in case to keep it
        break;
      }
      if (reord_idx == s) {
	reord_idx++;
	continue;
      }
      std::swap(matches_jk[reord_idx++], matches_jk[s]);
    }
    printf("IP: %d out of %d matches were reordered by model from %d to %d\n", reord_idx, matches_jk.size(), best_i, img_j);
    //if (reord_idx >= matches_jk.size()) {
     // break; //this part should be rewritten that prioritize points which we often have as inliers for other models
               //maybe just implement as a simple counter and order them in this way
    //}
    //break; //right now we will do only one iteratin for the testing purpose
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

void InlierPassing::store_matches_qual(image_t img_i, image_t img_j,
		                       std::vector<float> &matches_qual){
  m_qual[{img_i, img_j}] = matches_qual;
}

void InlierPassing::sort_matches_by_qual(const image_t img_i, const  image_t img_j,
					 FeatureMatches &matches){
  std::vector<std::pair<FeatureMatch, float>> merged(matches.size());
  std::vector<float> &q = m_qual[{img_i, img_j}];   
  for (size_t i = 0; i < merged.size(); ++i){
    merged[i].first = matches[i];
    merged[i].second = q[i];
  }
  //auto sortLambda = [](const auto &a, const auto &b) -> bool
  auto sortLambda = [](const std::pair<FeatureMatch, float> &a, const std::pair<FeatureMatch, float> &b) -> bool
                      { return a.second < b.second; };
  std::sort(merged.begin(), merged.end(), sortLambda);
  for (size_t i = 0; i < matches.size(); ++i){
    matches[i] = merged[i].first;
  }
  m_qual.erase({img_i, img_j});//clean memory, we don't need it anymore
}

std::unordered_map<image_t, std::vector<size_t>> 
const InlierPassing::get_inliers_passed(image_t img_j, image_t img_k) {
  std::unordered_map<image_t, std::vector<size_t>> res;
  if (connected_by.size() >= img_j + 1) {
    for (auto img_i : connected_by[img_j]) {
      //CHECK(img_i < img_j);   //not implemented for img_i > img_j
      //CHECK(img_j < img_k);   //not implemented for img_j > img_k
      //have to work without those checks
      //CHECK(pci.find[{img_i, img_j, img_k}] != pci.end());
      if (pci.find({img_i, img_j, img_k}) != pci.end()){      
        res[img_i] = pci[{img_i, img_j, img_k}];
      } else {
        //res[img_i] = std::vector<size_t>(0);	      
      }
    }
  }
  return res;
}

//bool InlierPassing::is_geom_found(image_t img1, image_t img2) { 
//  return static_cast<bool>(pair_inliers.count({img1, img2})); 
//};




} // namespace colmap
