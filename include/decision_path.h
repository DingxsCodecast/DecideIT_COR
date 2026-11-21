#ifndef DECIDEIT_DECISION_PATH_H
#define DECIDEIT_DECISION_PATH_H
#include "dblp.h"

namespace decision_path {
std::string sha256(const std::string& bytes);
std::string file_hash(const std::string& path);
std::string number(double x);
std::string quote(const std::string& x);
std::string array(const double* x, int n);
std::string rows_hash(const raw_data* p);
std::string bounds_hash(const raw_data* p);
double block_violation(const gsl_matrix* a, const gsl_vector* b,
                       const gsl_vector* lo, const gsl_vector* hi,
                       const double* point);
double point_violation(const raw_data* p, const double* probability,
                       const double* value, int first_dimension);
struct FocalInfo {
    int projection_lps = 0;
    double preparation_s = 0;
    double candidate_violation = 0;
    double final_violation = 0;
    bool loaded = false;
    std::string hash;
};
// Original Ax/Ay and bounds only; never relax statements to fit the focus.
FocalInfo feasible_focal(const raw_data* p, int first, int second,
                         double* probability, double* value,
                         const std::string& input_hash,
                         const std::string& load_path,
                         const std::string& save_path);
void init_log(const std::string& path);
void stage(const std::string& name, int level = -1);
void log(const std::string& json_fields);
void log_error(const std::string& message);
void direction(int level, const char* sense, const raw_data* model,
               int first_dimension, const double* probability,
               const double* value, const DecideItDirectionStats& stats,
               int aggregate_status, bool cplex, double bound,
               int cplex_status);
}
#endif
