#pragma once
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
inline void rlt_audit_snapshot(Highs &h) {
    const char* dir=std::getenv("DECIDEIT_RLT_SNAPSHOT_DIR");
    if(!dir)return;
    static unsigned count=0;
    if(count>=32)return;
    const std::string base=std::string(dir)+"/rlt_"+std::to_string(++count);
    highs_require(h.writeModel(base+".mps"),"audit model");
    highs_require(h.writeOptions(base+".options"),"audit options");
    if(h.getBasis().valid)highs_require(h.writeBasis(base+".basis"),"audit basis");
}
inline void rlt_audit_result(Highs &h,double seconds,bool fallback) {
    if(!std::getenv("DECIDEIT_RLT_DIAGNOSTICS"))return;
    const auto &b=h.getBasis();const auto &s=h.getSolution();const auto &lp=h.getLp();
    long long basic=0,near8=0,near7=0,near6=0;
    auto add=[&](double v,double lo,double hi){
        ++basic;double d=kHighsInf;
        if(std::isfinite(lo) && std::fabs(lo)<kHighsInf/2)d=std::fabs(v-lo);
        if(std::isfinite(hi) && std::fabs(hi)<kHighsInf/2)d=std::fmin(d,std::fabs(v-hi));
        if(d<=1e-8)++near8;if(d<=1e-7)++near7;if(d<=1e-6)++near6;
    };
    if(b.valid && s.value_valid){
        for(size_t j=0;j<b.col_status.size();++j)
            if(b.col_status[j]==HighsBasisStatus::kBasic)add(s.col_value[j],lp.col_lower_[j],lp.col_upper_[j]);
        for(size_t j=0;j<b.row_status.size();++j)
            if(b.row_status[j]==HighsBasisStatus::kBasic)add(s.row_value[j],lp.row_lower_[j],lp.row_upper_[j]);
    }
    const auto &i=h.getInfo();
    fprintf(stderr,"RLT_AUDIT {\"rows\":%d,\"columns\":%d,\"basis_valid\":%s,\"basic_variables\":%lld,\"basic_at_bound_1e8\":%lld,\"basic_at_bound_1e7\":%lld,\"basic_at_bound_1e6\":%lld,\"iterations\":%lld,\"max_primal_infeasibility\":%.17g,\"max_dual_infeasibility\":%.17g,\"solve_s\":%.9g,\"cold_fallback\":%s}\n",int(lp.num_row_),int(lp.num_col_),b.valid?"true":"false",basic,near8,near7,near6,(long long)i.simplex_iteration_count,i.max_primal_infeasibility,i.max_dual_infeasibility,seconds,fallback?"true":"false");
}
