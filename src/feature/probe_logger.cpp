//Ilia Shipachev, CVUT-FEL, CMP lab, 2020

#include "probe_logger.h"
#include <iostream>
#include <fstream>
#include "optim/ransac.h"
#include "optim/support_measurement.h"

namespace colmap {

  ProbeLogger::ProbeLogger(const std::string &filename, const bool isOneLiner_) :
    isOneLiner(isOneLiner_) {

    ostream.open(filename, std::ios::out);
    //ostream << "{";

    current_tab = isOneLiner ? " " : "\n";

    //ostream << inner_dict_start();
  }

  void ProbeLogger::write_header(const std::string &version, 
                                 const std::string &comment) {
    inner_dict_start();
    ostream << tab_kv_string("json_version", version);
    ostream << tab_kv_string("json_comment", comment);
  }

  void ProbeLogger::init_write_conf() {
    ostream << tab_key_string("conf");
    ostream << inner_dict_start();
  }

  void ProbeLogger::write_algo_conf(const std::unordered_map<std::string, int> &algo_conf) {
    ostream << tab_key_string("algo");
    ostream << inner_dict_start();
    for (auto key_val : algo_conf) {
      ostream << tab_kv_string(key_val.first, std::to_string(key_val.second));
    }
    ostream << inner_dict_end();
  }

  void ProbeLogger::write_sequential_matching_conf(const SequentialMatchingOptions &options) { 
    ostream << tab_key_string("sequential_matching");
    ostream << inner_dict_start();

    ostream << tab_kv_string("overlap", std::to_string(options.overlap));
  }

  void ProbeLogger::write_tvg_matching_conf(const SiftMatchingOptions& options) {
    ostream << tab_key_string("tvg_matching");
    ostream << inner_dict_start();

    ostream << tab_kv_string("num_threads", options.num_threads);
    ostream << tab_kv_string("use_gpu", options.use_gpu);
    ostream << tab_kv_string("max_ratio", options.max_ratio);
    ostream << tab_kv_string("max_distance", options.max_distance);
    ostream << tab_kv_string("max_error", options.max_error);
    ostream << tab_kv_string("confidence", options.confidence);
    ostream << tab_kv_string("min_num_trials", options.min_num_trials);
    ostream << tab_kv_string("max_num_trials", options.max_num_trials);
    ostream << tab_kv_string("min_inlier_ratio", options.min_inlier_ratio);
    ostream << tab_kv_string("guided_matching", options.guided_matching);
    ostream << tab_kv_string("inlier_passing", options.inlier_passing);
    ostream << inner_dict_end();
  }

  void ProbeLogger::finale_write_conf() {
    ostream << inner_dict_end();
  }

  void ProbeLogger::init_write_tvgs() {
    ostream << tab_key_string("tvgs");
    ostream << inner_arr_start();
  }

  void ProbeLogger::write_summary_tvg(image_t img1, image_t img2, int matches_num,
                                      int inl_num, int config, double time) {
    ostream << inner_dict_start();
    ostream << tab_kv_string("img_pair", img1, img2);
    ostream << tab_kv_string("matches_num",  matches_num);
    ostream << tab_kv_string("inl_num", inl_num);
    ostream << tab_kv_string("config", config);
    ostream << tab_kv_string("time", time);
  }

  template <typename Estimator, 
            typename Sampler>
  void ProbeLogger::write_model_report(typename const RANSAC<Estimator, 
                                                       InlierSupportMeasurer, 
                                                       Sampler>::Report &report, 
                                 const std::string model_type,
                                 double time) {
    ostream << tab_key_string(model_type);  //F or H or E
    ostream << inner_dict_start();
    ostream << tab_kv_string("num_trials", report.num_trials);  //F or H or E
    ostream << tab_kv_string("num_inliers", report.support.num_inliers);
    ostream << tab_kv_string("time", time);
    ostream << inner_dict_end();
  }

  void ProbeLogger::write_inl_passed_stat(image_t img_id,
                                          const std::unordered_map<image_t, size_t> &inl_passed) {
    ostream << tab_key_string("inl_passing");
    ostream << inner_arr_start();
    for (auto key_val : inl_passed) {
      ostream << inner_dict_start();
      ostream << tab_kv_string("from", static_cast<int>(key_val.first));
      ostream << tab_kv_string("pci", static_cast<int>(key_val.second));
      ostream << inner_dict_end();
    }
    ostream << inner_arr_end();
  }

  void ProbeLogger::finale_write_tvgs() {
    ostream << inner_dict_end();
  }

  void ProbeLogger::finile_write_all() {
    ostream << inner_dict_end();
  }

  ProbeLogger::~ProbeLogger() {
    ostream.close();
  }

  //PRIVATE FUNCTIONS

  std::string ProbeLogger::inner_dict_start() {
    std::string res;
    res = current_tab + std::string("{");
    if (!isOneLiner) {
      current_tab.append("\t");
    }
    return res;
  }

  std::string ProbeLogger::inner_arr_start() {
    std::string res;
    res = current_tab + std::string("[");
    if (!isOneLiner) {
      current_tab.append("\t");
    }
    return res;
  }

  std::string ProbeLogger::inner_dict_end() {
    if (!isOneLiner) {
      current_tab.pop_back();
    }
    return current_tab + std::string("}") + std::string(",");
  }

  std::string ProbeLogger::inner_arr_end() {
    if (!isOneLiner) {
      current_tab.pop_back();
    }
    return current_tab + std::string("]") + std::string(",");
  }

  std::string ProbeLogger::tab_kv_string(const std::string &key, const std::string &val) {
    return current_tab + key + std::string(": ") + val + std::string(",");
  }

  std::string ProbeLogger::tab_kv_string(const std::string &key, int val) {
    return current_tab + key + std::string(": ") + std::to_string(val) + std::string(",");
  }

  std::string ProbeLogger::tab_kv_string(const std::string &key, double val) {
    return current_tab + key + std::string(": ") + std::to_string(val) + std::string(",");
  }

  std::string ProbeLogger::tab_kv_string(const std::string &key, int val1, int val2) {
    return (  current_tab + key + std::string(": ") 
            + std::string("[") + std::to_string(val1) + std::string(", ")
            + std::to_string(val2) + std::string("]")
            + std::string(",")
           );
  }

  std::string ProbeLogger::tab_key_string(const std::string &key) {
    return current_tab + key + std::string(":");
  }

}