#ifndef DECIDEIT_RESOURCE_SCOPE_H
#define DECIDEIT_RESOURCE_SCOPE_H
#include <unordered_map>
#include <new>

// Retain the existing GSL-based geometry while making early returns and
// exceptions release per-query allocations. Explicit frees unregister objects.
// Views are never registered; only owning allocation functions pass through here.
namespace managed {
struct Scope {
    using Deleter=void(*)(void*);
    static thread_local Scope* current;
    Scope* previous;
    std::unordered_map<void*,Deleter> live;
    Scope():previous(current){current=this;}
    Scope(const Scope&)=delete;
    Scope& operator=(const Scope&)=delete;
    ~Scope(){while(!live.empty()){auto i=live.begin();auto p=i->first;auto d=i->second;live.erase(i);d(p);}current=previous;}
    static void forget(void* p){for(auto s=current;s;s=s->previous)if(s->live.erase(p))break;}
};
template<class T,void(*Release)(T*)> T* own(T* p){
    if(!p)throw std::bad_alloc();
    if(Scope::current){try{Scope::current->live.emplace(p,[](void* v){Release(static_cast<T*>(v));});}catch(...){Release(p);throw;}}
    return p;
}
template<class T,void(*Release)(T*)> void drop(T* p){if(p){Scope::forget(p);Release(p);}}
#define DECIDEIT_GSL_OWNER(name,type) \
inline type* name##_alloc(size_t n){return own<type,::gsl_##name##_free>(::gsl_##name##_alloc(n));} \
inline type* name##_calloc(size_t n){return own<type,::gsl_##name##_free>(::gsl_##name##_calloc(n));} \
inline void name##_free(type* p){drop<type,::gsl_##name##_free>(p);}
DECIDEIT_GSL_OWNER(vector,gsl_vector)
DECIDEIT_GSL_OWNER(vector_int,gsl_vector_int)
DECIDEIT_GSL_OWNER(permutation,gsl_permutation)
#undef DECIDEIT_GSL_OWNER
#define DECIDEIT_GSL_MATRIX(name,type) \
inline type* name##_alloc(size_t m,size_t n){return own<type,::gsl_##name##_free>(::gsl_##name##_alloc(m,n));} \
inline type* name##_calloc(size_t m,size_t n){return own<type,::gsl_##name##_free>(::gsl_##name##_calloc(m,n));} \
inline void name##_free(type* p){drop<type,::gsl_##name##_free>(p);}
DECIDEIT_GSL_MATRIX(matrix,gsl_matrix)
DECIDEIT_GSL_MATRIX(matrix_int,gsl_matrix_int)
#undef DECIDEIT_GSL_MATRIX
template<class T> void delete_object(T* p){delete p;}
template<class T> T* create(){return own<T,delete_object<T>>(new T);}
template<class T> void destroy(T* p){drop<T,delete_object<T>>(p);}
}
#endif
