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
    //0tab, begin of the file
    ProbeLogger(const std::string &filename, const bool isOneLine);
    //void init_write_all();

    //1tab -- header -- BEGIN
    void write_head_open(const std::string &version, const std::string &comment); //writes header of the version and some comment about the application

    //2tab -- configs -- BEGIN
    void write_conf_open();
    void write_conf_algo(const std::unordered_map<std::string, int> &algo_conf);
    void write_sequential_matching_conf(const SequentialMatchingOptions &options);
    void write_tvg_matching_conf(const SiftMatchingOptions& options);
    void write_inl_passed_stat(image_t img_id, const std::unordered_map<image_t, size_t> &inl_passed);
    void write_conf_close();
    //2tab -- configs -- END

    //2tab -- tvgs -- BEGIN
    void write_tvgs_open();   //two view geometries
    void write_summary_tvg_open(image_t img1, image_t img2, int matches_num,
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
    void write_tvgs_close();
    //2tab -- tvgs -- END
    
    void write_head_close();
    //1tab -- header -- END

    //0tab -- end of the file
    ~ProbeLogger();
    
  private:
    bool isOneLiner = false;
    std::string current_tab = "";
  
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

  };

}  // namespace colmap
