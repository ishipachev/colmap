//Ilia Shipachev, CVUT-FEL, CMP lab, 2020

#pragma once
#include "util/types.h"
#include "feature/types.h"
#include <unordered_map>
#include <vector>
#include <string>
#include "feature/matching.h"
#include "feature/sift.h"
#include <fstream>
#include <iostream>
#include "optim/ransac.h"
#include "optim/support_measurement.h"
#include "feature/inlier_passing.h"

#define INL_PASSING_ON false

namespace colmap {

  class ProbeLogger {

  public:
    ProbeLogger(const std::string &filename, const bool isOneLine);
    void init_write_all();

    //1tab -- header -- BEGIN
    void write_header(const std::string &version, const std::string &comment); //writes header of the version and some comment about the application

    //2tab -- configs -- BEGIN
    void init_write_conf();
    void write_algo_conf(const std::unordered_map<std::string, int> &algo_conf);
    void write_sequential_matching_conf(const SequentialMatchingOptions &options);
    void write_tvg_matching_conf(const SiftMatchingOptions& options);
    void write_inl_passed_stat(image_t img_id, 
                               const std::unordered_map<image_t, size_t> &inl_passed);
    void finale_write_conf();
    //2tab -- configs -- END

    //2tab -- tvgs -- BEGIN
    void init_write_tvgs();   //two view geometries
    void write_summary_tvg(image_t img1, image_t img2, int matches_num,
                           int inl_num, int config, double time);

    template <typename Estimator, 
              typename Sampler>
    void write_model_report(typename const RANSAC<Estimator, 
                                              InlierSupportMeasurer,
                                              Sampler>::Report &report,
                        const std::string model_type,
                        double time);
    void write_inl_passed_stat(image_t img_id,
                               const std::unordered_map<image_t, size_t> &inl_passed);

    void finale_write_tvgs();
    //2tab -- tvgs -- END
    
    void finile_write_all();


    ~ProbeLogger();
    
  private:

    int current_tab_cnt = 0;
    std::string current_tab = "";
    
    bool isOneLiner = false;
    std::ofstream ostream;

    std::string inner_dict_start();
    std::string inner_dict_end();
    std::string inner_arr_start();
    std::string inner_arr_end();

    std::string ProbeLogger::tab_kv_string(const std::string &key, const std::string &val);
    std::string ProbeLogger::tab_kv_string(const std::string &key, int val);
    std::string ProbeLogger::tab_kv_string(const std::string &key, double val);
    std::string ProbeLogger::tab_kv_string(const std::string &key, int val1, int val2);

    std::string ProbeLogger::tab_key_string(const std::string &key);

    //std::unordered_map<std::pair<image_t, image_t>, FeatureMatches> pair_inliers; //inliers of model i-to-j
    std::unordered_map<std::pair<image_t, image_t>, std::vector<point2D_t>> pair_inliers; //j's inliers of model i-to-j

    std::vector<std::vector<image_t>> connected_by; //list of all models builded to i-s image
  };

}  // namespace colmap
