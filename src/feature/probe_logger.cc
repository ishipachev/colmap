//Ilia Shipachev, CVUT-FEL, CMP lab, 2020

#include "probe_logger.h"
#include <iostream>
#include <fstream>
#include <utility>
#include "optim/ransac.h"
#include "optim/support_measurement.h"
#include "estimators/essential_matrix.h"
#include "estimators/fundamental_matrix.h"
#include "estimators/homography_matrix.h"

namespace colmap {

  ProbeLogger probeLogger;    //bringing probeLogger to the global colmap namespace

  ProbeLogger::ProbeLogger() {
  }

  void ProbeLogger::init(const std::string &filename, const bool isOneLiner_) {
    isOneLiner = isOneLiner_;
    printf("ProbeLogger filename: %s \n", filename);
    ostream.open(filename);
    current_tab = isOneLiner ? " " : "\n";
  }

  void ProbeLogger::write_head_open(const std::string &version, 
                                    const std::string &comment) {
    ostream << inner_dict_start();
    ostream << tab_kv_string("json_version", version);
    ostream << tab_kv_string("json_comment", comment);
  }

  void ProbeLogger::write_conf_open() {
    ostream << tab_key_string("conf");
    ostream << inner_dict_start();
  }

  void ProbeLogger::write_conf_algo(const std::unordered_map<std::string, int> &algo_conf) {
    ostream << tab_key_string("algo");
    ostream << inner_dict_start();
    for (auto key_val : algo_conf) {
      ostream << tab_kv_string(key_val.first, key_val.second);
    }
    ostream << inner_dict_end();
  }

  void ProbeLogger::write_conf_sequential_matching(const struct SequentialMatchingOptions &options) { 
    ostream << tab_key_string("sequential_matching");

    ostream << inner_dict_start();
    ostream << tab_kv_string("overlap", std::to_string(options.overlap));
    ostream << inner_dict_end();
  }

  void ProbeLogger::write_conf_tvg_matching(const struct SiftMatchingOptions& options) {
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

  void ProbeLogger::write_conf_close() {
    ostream << inner_dict_end();
  }

  void ProbeLogger::write_tvgs_open() {
    ostream << tab_key_string("tvgs");
    ostream << inner_arr_start();
  }

  void ProbeLogger::write_tvg_open(image_t img1, image_t img2, int matches_num) {
    ostream << inner_dict_start();
    ostream << tab_kv_string("img_pair", img1, img2);
    ostream << tab_kv_string("matches_num",  matches_num);
  }

  void ProbeLogger::store_matches_qual(image_t img1, image_t img2, 
                                       std::vector<float> &matches_qual) {
    m_qual[{img1, img2}] = std::move(matches_qual);
  }

  void ProbeLogger::write_stored_matches_qual(image_t img1, image_t img2) {
    //ostream << tab_kvs_string("matches_quality", m_qual[{img1, img2}]);
    //to decrease the size of logs have it commented out
    ostream << tab_kv_string("matches_quality:","[]");
    m_qual.erase({img1, img2});  //we don't need this array after printing
  }

  void ProbeLogger::write_ransac_open() {
    ostream << tab_key_string("inl_cnt") << std::string("[");
  }

  void ProbeLogger::write_ransac_inlier_cnt(size_t cnt) {
    ostream << std::to_string(cnt) << std::string(", ");
  }

  void ProbeLogger::write_ransac_close() {
    ostream << "]";
  }

  void ProbeLogger::write_model_report_open(const std::string model_type) {
    ostream << tab_key_string(model_type);  //F or H or E
    ostream << inner_dict_start();
  }

  template <typename Estimator, 
            typename SupportMeasurer>
  void ProbeLogger::write_model_report(const typename RANSAC<Estimator, 
                                                             SupportMeasurer>::Report &report, 
                                       //const std::string model_type,
                                       double time) {
    //ostream << tab_key_string(model_type);  //F or H or E
    //ostream << inner_dict_start();
    ostream << tab_kv_string("num_trials", report.num_trials);  //F or H or E
    ostream << tab_kv_string("num_inliers", report.support.num_inliers);
    ostream << tab_kv_string("time", time);
    //ostream << inner_dict_end();
  }

  void ProbeLogger::write_model_report_close() {
    ostream << inner_dict_end();
  }

  void ProbeLogger::write_inl_passed_stat(const std::unordered_map<image_t, size_t> &inl_passed) {
    ostream << tab_key_string("inl_passed");
    ostream << inner_arr_start();
    for (auto key_val : inl_passed) {
      ostream << inner_dict_start();
      ostream << tab_kv_string("from", static_cast<size_t>(key_val.first));
      ostream << tab_kv_string("pci", static_cast<size_t>(key_val.second));
      ostream << inner_dict_end();
    }
    ostream << inner_arr_end();
  }

  void ProbeLogger::write_inl_passed_stat(const std::unordered_map<image_t, std::vector<size_t>> &inl_passed) {
    ostream << tab_key_string("inl_passed");
    ostream << inner_arr_start();
    for (auto key_val : inl_passed) {
      ostream << inner_dict_start();
      ostream << tab_kv_string("from", static_cast<size_t>(key_val.first));
      ostream << tab_kv_string("pci", static_cast<size_t>(key_val.second.size()));
      //to decrease the size of logs commented out array logs
      //ostream << tab_kvs_string("inds", key_val.second);
      ostream << tab_kv_string("inds:", "[]");
      ostream << inner_dict_end();
    }
    ostream << inner_arr_end();
  }

  void ProbeLogger::write_inliers(const std::vector<char>* const inlier_mask) {
    std::vector<size_t> inds(inlier_mask->size());
    int cnt=0;
    for (int i=0; i < inlier_mask->size(); ++i) {
      if ((*inlier_mask)[i]) {
        inds[cnt++] = i;
      }
    }
    inds.resize(cnt);
    //ostream << tab_kvs_string("inliers", inds);
    //to decrease the size of logs commented out array logs
    ostream << tab_kv_string("inliers:", "[]");
  }

  void ProbeLogger::write_tvg_close(size_t inl_num, int config, double time) {
    ostream << tab_kv_string("inl_num", static_cast<int>(inl_num));
    ostream << tab_kv_string("config", config);
    ostream << tab_kv_string("time", time);
    ostream << inner_dict_end();
  }

  void ProbeLogger::write_tvgs_close() {
    ostream << inner_arr_end();
  }

  void ProbeLogger::write_head_close() {
    std::string last_line;
    last_line = inner_dict_end();
    last_line.pop_back(); //getting rid of trailing comma
    ostream << last_line;
  }

  void ProbeLogger::deinit() {
    //ostream << inner_dict_end();
    ostream.close();
  }

  ProbeLogger::~ProbeLogger() {
    //inner_dict_end();
    //ostream.close();
  }

  //PRIVATE FUNCTIONS

  std::string ProbeLogger::inner_dict_start() {
    std::string res;
    res = current_tab + std::string("{");
    if (!isOneLiner) {
      current_tab.append(tab_str);
    }
    return res;
  }

  std::string ProbeLogger::inner_arr_start() {
    std::string res;
    res = current_tab + std::string("[");
    if (!isOneLiner) {
      current_tab.append(tab_str);
    }
    return res;
  }

  std::string ProbeLogger::inner_dict_end() {
    if (!isOneLiner) {
      //printf("Debug print:\n");
      //printf("\t current_tab=%s, size=%d\n;\t tab_str=%s, size=%d\n",
		//current_tab, current_tab.size(), tab_str, tab_str.size());
      current_tab.erase(current_tab.end() - tab_str.size());
      //current_tab.pop_back(); //was done for '\t' symbol
    }
    return current_tab + std::string("}") + std::string(",");
  }

  std::string ProbeLogger::inner_arr_end() {
    if (!isOneLiner) {
      current_tab.erase(current_tab.size() - tab_str.size() - 1, tab_str.size());
      //current_tab.pop_back(); //was done for '\t' symbol
    }
    return current_tab + std::string("]") + std::string(",");
  }

  std::string ProbeLogger::tab_kv_string(const std::string &key, const std::string &val) {
    return current_tab + std::string("\"") + key + std::string("\": ") 
         + std::string("\"") + val + std::string("\",");
  }

  std::string ProbeLogger::tab_kv_string(const std::string &key, int val) {
    return current_tab + std::string("\"") + key + std::string("\": ")
         + std::to_string(val) + std::string(",");
  }

  std::string ProbeLogger::tab_kv_string(const std::string &key, size_t val) {
    return current_tab + std::string("\"") + key + std::string("\": ") 
         + std::to_string(val) + std::string(",");
  }

  std::string ProbeLogger::tab_kv_string(const std::string &key, double val) {
    return current_tab + std::string("\"") + key + std::string("\": ") 
         + std::to_string(val) + std::string(",");
  }

  std::string ProbeLogger::tab_kv_string(const std::string &key, int val1, int val2) {
    return (  current_tab + std::string("\"") + key + std::string("\": ") 
            + std::string("[") + std::to_string(val1) + std::string(", ")
            + std::to_string(val2) + std::string("],")
           );
  }

  std::string ProbeLogger::tab_kvs_string(const std::string &key, const std::vector<float> &vals) {
    std::string res;
    res = current_tab + std::string("\"") + key + std::string("\": [");
    char s_buf[8];
    for (auto v : vals) {
      std::sprintf(s_buf, "%4.3f", v);
      res += std::string(s_buf) + std::string(", ");
    }
    res += std::string("],");
    return res;
  }

  std::string ProbeLogger::tab_kvs_string(const std::string &key, const std::vector<size_t> &vals) {
    std::string res;
    res = current_tab + std::string("\"") + key + std::string("\": [");
    for (auto v : vals) {
      res += std::to_string(v) + std::string(", ");
    }
    res += std::string("],");
    return res;
  }

  std::string ProbeLogger::tab_key_string(const std::string &key) {
    return current_tab + std::string("\"") + key + std::string("\":");
  }

  template void ProbeLogger::write_model_report<
    EssentialMatrixFivePointEstimator,
    InlierSupportMeasurer>(const typename RANSAC<EssentialMatrixFivePointEstimator,
                           InlierSupportMeasurer>::Report &report,
                           //const std::string model_type,
                           double time);

  template void ProbeLogger::write_model_report<
    FundamentalMatrixSevenPointEstimator,
    InlierSupportMeasurer>(const typename RANSAC<FundamentalMatrixSevenPointEstimator,
                           InlierSupportMeasurer>::Report &report,
                           //const std::string model_type,
                           double time);

  template void ProbeLogger::write_model_report<
    HomographyMatrixEstimator,
    InlierSupportMeasurer>(const typename RANSAC<HomographyMatrixEstimator,
                           InlierSupportMeasurer>::Report &report,
                           //const std::string model_type,
                           double time);
}


