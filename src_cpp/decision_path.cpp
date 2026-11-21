#include "decision_path.h"
#include <array>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cerrno>
#include <fcntl.h>
#include <unistd.h>

namespace decision_path {
namespace {
FILE* diagnostics = nullptr;
std::string current_stage = "input";
int current_level = -1;
uint32_t rotate(uint32_t x, int n) { return (x >> n) | (x << (32-n)); }
void exclusive_write(const std::string& path, const std::string& bytes) {
    const int fd = ::open(path.c_str(), O_WRONLY|O_CREAT|O_EXCL, 0644);
    if (fd < 0) throw std::runtime_error("Cannot create new artifact: " + path);
    size_t done = 0;
    while (done < bytes.size()) {
        const ssize_t n = ::write(fd, bytes.data()+done, bytes.size()-done);
        if (n < 0 && errno == EINTR) continue;
        if (n <= 0) { ::close(fd); throw std::runtime_error("Artifact write failed: "+path); }
        done += static_cast<size_t>(n);
    }
    ::fsync(fd);
    ::close(fd);
}
// L1 projection is needed only for an infeasible candidate. The two
// probability normalizations are already present as pairs of rows in Ay.
void project(const gsl_matrix* a, const gsl_vector* b,
             const gsl_vector* lo, const gsl_vector* hi, double* point) {
    const int n = static_cast<int>(lo->size);
    HighsLp lp;
    lp.num_col_ = 2*n;
    lp.num_row_ = static_cast<HighsInt>(a->size1)+2*n;
    lp.col_cost_.assign(2*n, 0.0);
    lp.col_lower_.assign(2*n, 0.0);
    lp.col_upper_.assign(2*n, kHighsInf);
    for (int j=0; j<n; ++j) {
        lp.col_lower_[j]=gsl_vector_get(lo,j);
        lp.col_upper_[j]=gsl_vector_get(hi,j);
        lp.col_cost_[n+j]=1.0;
    }
    lp.row_lower_.assign(lp.num_row_, -kHighsInf);
    lp.row_upper_.resize(lp.num_row_);
    lp.a_matrix_.format_=MatrixFormat::kRowwise;
    lp.a_matrix_.start_.assign(1,0);
    for (size_t r=0; r<a->size1; ++r) {
        lp.row_upper_[r]=gsl_vector_get(b,r);
        for (int j=0; j<n; ++j) {
            const double c=gsl_matrix_get(a,r,j);
            if (c!=0) {lp.a_matrix_.index_.push_back(j); lp.a_matrix_.value_.push_back(c);}
        }
        lp.a_matrix_.start_.push_back(lp.a_matrix_.value_.size());
    }
    for (int j=0; j<n; ++j) for (int s=-1; s<=1; s+=2) {
        const int r=static_cast<int>(a->size1)+2*j+(s==1);
        lp.row_upper_[r]=s*point[j+1];
        lp.a_matrix_.index_.push_back(j); lp.a_matrix_.value_.push_back(s);
        lp.a_matrix_.index_.push_back(n+j); lp.a_matrix_.value_.push_back(-1.0);
        lp.a_matrix_.start_.push_back(lp.a_matrix_.value_.size());
    }
    std::lock_guard<std::recursive_mutex> lock(highs_api_mutex());
    Highs highs;
    configure_highs(highs);
    highs_require(highs.passModel(lp), "focal L1 model");
    highs_require(highs.run(), "focal L1 solve");
    if (highs.getModelStatus()!=HighsModelStatus::kOptimal ||
        !highs_solution_valid(highs))
        throw std::runtime_error("Original block infeasible or focal projection uncertified: status="+
            std::to_string(static_cast<int>(highs.getModelStatus()))+", primal="+
            number(highs.getInfo().max_primal_infeasibility)+", dual="+
            number(highs.getInfo().max_dual_infeasibility));
    const auto& sol=highs.getSolution();
    for (int j=0; j<n; ++j) point[j+1]=sol.col_value[j];
}
std::string focal_payload(const raw_data* p, int first, int second,
                          const double* probability, const double* value,
                          const std::string& input_hash) {
    std::ostringstream s;
    s.imbue(std::locale::classic());
    s << "DECIDEIT_FEASIBLE_FOCAL_V1\n" << input_hash << '\n'
      << rows_hash(p) << '\n' << bounds_hash(p) << '\n'
      << first << ' ' << second << '\n' << std::setprecision(17);
    for(int j=1;j<=first+second;++j) s<<probability[j]<<' ';
    s<<'\n';
    for(int j=1;j<=first+second;++j) s<<value[j]<<' ';
    s<<'\n';
    return s.str();
}
}
// SHA-256 of exact bytes. Model hashes use a canonical, classic-locale,
// round-trip decimal encoding; input hashes use the actual DDT bytes.
std::string sha256(const std::string& bytes) {
    static const uint32_t k[64]={
      0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
      0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
      0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
      0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
      0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
      0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
      0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
      0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
    uint32_t h[8]={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    std::vector<unsigned char> data(bytes.begin(),bytes.end());
    const uint64_t bits=static_cast<uint64_t>(data.size())*8;
    data.push_back(0x80);
    while(data.size()%64!=56) data.push_back(0);
    for(int j=7;j>=0;--j) data.push_back(static_cast<unsigned char>(bits>>(8*j)));
    for(size_t pos=0;pos<data.size();pos+=64) {
        uint32_t w[64];
        for(int j=0;j<16;++j) {
            w[j]=0;
            for(int z=0;z<4;++z) w[j]=(w[j]<<8)|data[pos+4*j+z];
        }
        for(int j=16;j<64;++j) {
            const uint32_t x=w[j-15], y=w[j-2];
            w[j]=w[j-16]+(rotate(x,7)^rotate(x,18)^(x>>3))+w[j-7]+(rotate(y,17)^rotate(y,19)^(y>>10));
        }
        uint32_t a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],v=h[7];
        for(int j=0;j<64;++j) {
            const uint32_t t1=v+(rotate(e,6)^rotate(e,11)^rotate(e,25))+((e&f)^(~e&g))+k[j]+w[j];
            const uint32_t t2=(rotate(a,2)^rotate(a,13)^rotate(a,22))+((a&b)^(a&c)^(b&c));
            v=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
        }
        h[0]+=a;h[1]+=b;h[2]+=c;h[3]+=d;h[4]+=e;h[5]+=f;h[6]+=g;h[7]+=v;
    }
    std::ostringstream s;
    s<<std::hex<<std::setfill('0');
    for(int j=0;j<8;++j) s<<std::setw(8)<<h[j];
    return s.str();
}
std::string file_hash(const std::string& path) {
    std::ifstream f(path.c_str(),std::ios::binary);
    if(!f) throw std::runtime_error("Cannot hash input: "+path);
    std::ostringstream s; s<<f.rdbuf(); return sha256(s.str());
}
std::string number(double x) {
    if(!std::isfinite(x)) return "null";
    std::ostringstream s; s.imbue(std::locale::classic()); s<<std::setprecision(17)<<x; return s.str();
}
std::string quote(const std::string& x) {
    std::ostringstream s; s<<'"';
    for(unsigned char c:x) {
        if(c=='"'||c=='\\') s<<'\\'<<c;
        else if(c<32) s<<"\\u00"<<std::hex<<std::setw(2)<<std::setfill('0')<<static_cast<int>(c)<<std::dec;
        else s<<c;
    }
    s<<'"'; return s.str();
}
std::string array(const double* x,int n) {
    std::string s="[";
    for(int j=1;j<=n;++j) {if(j>1)s+=",";s+=number(x[j]);}
    return s+"]";
}
std::string rows_hash(const raw_data* p) {
    std::ostringstream s; s.imbue(std::locale::classic()); s<<std::setprecision(17);
    const gsl_matrix* aa[2]={p->Ay,p->Ax}; const gsl_vector* bb[2]={p->by,p->bx};
    for(int k=0;k<2;++k) {
        s<<aa[k]->size1<<' '<<aa[k]->size2<<'\n';
        for(size_t i=0;i<aa[k]->size1;++i) {
            for(size_t j=0;j<aa[k]->size2;++j) s<<gsl_matrix_get(aa[k],i,j)<<' ';
            s<<gsl_vector_get(bb[k],i)<<'\n';
        }
    }
    return sha256(s.str());
}
std::string bounds_hash(const raw_data* p) {
    std::ostringstream s; s.imbue(std::locale::classic()); s<<std::setprecision(17);
    const gsl_vector* v[4]={p->yl,p->yu,p->xl,p->xu};
    for(int k=0;k<4;++k) {for(size_t j=0;j<v[k]->size;++j)s<<gsl_vector_get(v[k],j)<<' ';s<<'\n';}
    return sha256(s.str());
}
double block_violation(const gsl_matrix* a,const gsl_vector* b,
                       const gsl_vector* lo,const gsl_vector* hi,const double* point) {
    double result=0;
    for(size_t j=0;j<lo->size;++j) {
        if(!std::isfinite(point[j+1])) return kHighsInf;
        result=std::fmax(result,gsl_vector_get(lo,j)-point[j+1]);
        result=std::fmax(result,point[j+1]-gsl_vector_get(hi,j));
    }
    for(size_t r=0;r<a->size1;++r) {
        double lhs=0;
        for(size_t j=0;j<a->size2;++j) lhs+=gsl_matrix_get(a,r,j)*point[j+1];
        result=std::fmax(result,lhs-gsl_vector_get(b,r));
    }
    return result;
}
double point_violation(const raw_data* p,const double* prob,const double* val,int first) {
    double result=std::fmax(block_violation(p->Ay,p->by,p->yl,p->yu,prob),
                            block_violation(p->Ax,p->bx,p->xl,p->xu,val));
    double a=0,b=0;
    for(size_t j=1;j<=p->yl->size;++j) (j<=static_cast<size_t>(first)?a:b)+=prob[j];
    return std::fmax(result,std::fmax(std::fabs(a-1),std::fabs(b-1)));
}
FocalInfo feasible_focal(const raw_data* p,int first,int second,double* prob,double* val,
                         const std::string& input_hash,const std::string& load_path,const std::string& save_path) {
    const auto begin=std::chrono::steady_clock::now();
    FocalInfo info;
    info.candidate_violation=point_violation(p,prob,val,first);
    if(!load_path.empty()) {
        std::ifstream f(load_path.c_str());
        std::string magic,input,rows,bounds,hash; int a=0,b=0;
        if(!(f>>magic>>input>>rows>>bounds>>a>>b) || magic!="DECIDEIT_FEASIBLE_FOCAL_V1" ||
           input!=input_hash || rows!=rows_hash(p) || bounds!=bounds_hash(p) || a!=first || b!=second)
            throw std::runtime_error("Focal artifact/input/model mismatch: "+load_path);
        for(int j=1;j<=first+second;++j) f>>prob[j];
        for(int j=1;j<=first+second;++j) f>>val[j];
        f>>hash>>info.projection_lps>>info.preparation_s;
        std::string extra;
        if(!f || f>>extra || hash!=sha256(focal_payload(p,first,second,prob,val,input_hash)) ||
           info.projection_lps<0 || info.projection_lps>2 || !std::isfinite(info.preparation_s) || info.preparation_s<0)
            throw std::runtime_error("Invalid focal artifact: "+load_path);
        info.loaded=true;
    } else {
        if(block_violation(p->Ay,p->by,p->yl,p->yu,prob)>1e-9) {
            project(p->Ay,p->by,p->yl,p->yu,prob); ++info.projection_lps;
        }
        if(block_violation(p->Ax,p->bx,p->xl,p->xu,val)>1e-9) {
            project(p->Ax,p->bx,p->xl,p->xu,val); ++info.projection_lps;
        }
        info.preparation_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    }
    info.final_violation=point_violation(p,prob,val,first);
    if(!std::isfinite(info.final_violation) || info.final_violation>1e-8)
        throw std::runtime_error("Focus is not feasible in ORIGINAL decision model");
    const std::string payload=focal_payload(p,first,second,prob,val,input_hash);
    info.hash=sha256(payload);
    if(!save_path.empty()) exclusive_write(save_path,payload+info.hash+"\n"+
                       std::to_string(info.projection_lps)+" "+number(info.preparation_s)+"\n");
    return info;
}
void init_log(const std::string& path) {
    setvbuf(stdout,nullptr,_IOLBF,0);
    if(path.empty()) return;
    const int fd=::open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL,0644);
    if(fd<0) throw std::runtime_error("Diagnostics file must be new: "+path);
    diagnostics=fdopen(fd,"w");
    if(!diagnostics) {::close(fd);throw std::runtime_error("Cannot open diagnostics stream");}
}
void stage(const std::string& name,int level) {current_stage=name;current_level=level;}
void log(const std::string& fields) {
    const std::string record="{\"schema\":\"decideit-diagnostics-v1\","+fields+"}";
    printf("DIAGNOSTIC %s\n",record.c_str());
    if(diagnostics) {fprintf(diagnostics,"%s\n",record.c_str());fflush(diagnostics);}
}
void log_error(const std::string& message) {
    log("\"event\":\"error\",\"stage\":"+quote(current_stage)+",\"level\":"+
        std::to_string(current_level)+",\"message\":"+quote(message));
}
void direction(int level,const char* sense,const raw_data* model,int first,
               const double* prob,const double* val,const DecideItDirectionStats& stats,
               int aggregate_status,bool cplex,double bound,int cplex_status) {
    const double violation=point_violation(model,prob,val,first);
    double obj=0;for(size_t j=1;j<=model->yl->size;++j)obj+=(j<=static_cast<size_t>(first)?1:-1)*prob[j]*val[j];
    const double reported=std::string(sense)=="min"?prob[0]:-prob[0];
    const bool attempted = cplex ? cplex_status!=0 :
        aggregate_status==OK || stats.certificate!="not_recorded" || stats.amc_rounds>0;
    log("\"event\":\"direction\",\"level\":"+std::to_string(level)+",\"sense\":"+quote(sense)+
        ",\"attempted\":"+(attempted?"true":"false")+
        ",\"aggregate_status\":"+std::to_string(aggregate_status)+",\"termination_flag\":"+
        (cplex?"null":std::to_string(stats.termination_flag))+
        ",\"cplex_status\":"+(cplex?std::to_string(cplex_status):"null")+
        ",\"incumbent\":"+number(reported)+",\"witness_objective\":"+number(obj)+
        ",\"witness_violation\":"+number(violation)+",\"objective_residual\":"+number(std::fabs(obj-reported))+
        ",\"bound\":"+number(bound)+",\"bound_scope\":"+quote(cplex?"original_level":"residual_after_valid_cuts")+
        ",\"certificate\":"+quote(cplex?"cplex_status":stats.certificate)+
        ",\"rlt_primal_infeasibility\":"+number(stats.rlt_primal_infeasibility)+
        ",\"rlt_dual_infeasibility\":"+number(stats.rlt_dual_infeasibility)+
        ",\"cuts\":"+std::to_string(stats.cuts)+",\"p\":"+array(prob,model->yl->size)+
        ",\"v\":"+array(val,model->xl->size));
}
}
