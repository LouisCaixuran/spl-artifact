// To use this from a port do the following:
//
// 1) Supply a cost function
// template <class G_t, class I_t>
// float instruction_cost(const assignment &a, unsigned short int i, const G_t &G, const I_t &I);
// Which can range from
// simple, e.g. cost 1 for each byte accessed in a register, cost 4 for each byte accessed in memory
// to
// quite involved, e.g. the number of bytes of code the code generator would generate.
//1)create a cost function 
// float instruction_cost(const i_assignment_ps &a)
//
// 2) Call
// create_cfg(), convert_cfg_to_spcfg() and generate_spcfg() in that order


//TODO: change the global_reg to map

#include <boost/version.hpp>
#if BOOST_VERSION == 106000
   #include <boost/type_traits/ice.hpp>
#endif

#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/container/flat_map.hpp>

#include "common.h"
#include "SDCCcfg.hpp"
#include <chrono>


extern "C"
{
#include "SDCCbtree.h"
}

// Integer constant upper bound on port->num_regs
#define MAX_NUM_REGS 3
//int duration_of_permutation=0;

static int getIndex(std::vector<short int> v, short int K) 
{ 
    auto it = find(v.begin(), v.end(), K); 
  
    // If element was found 
    if (it != v.end())  
    { 
      
        // calculating the index 
        // of K 
        int index = it - v.begin(); 
        return index;
    } 
    else { 
        // If the element is not 
        // present in the vector 
        return -1;
    } 
} 


static void convert_to_global(std::vector<short int> v,std::vector<var_t> variables, f &global, int n){ 
   global.resize(n,-3);
   int end=variables.size();
   for(int i=0;i<end;++i){
      global[variables[i]]=getIndex(v,variables[i]);
   }}


static void if_f_match(f f1,f f2, f &f3){
   int n=f1.size();
  
   f3.resize(n);
   for(int i=0;i<n;++i){
      if(f1[i]==-3||f1[i]==f2[i]){
         f3[i]=f2[i];
      }else if (f2[i]==-3){
         f3[i]=f1[i];
      }else{
         f3[0]=-2;
         return;
      
      }
   }
}

//get to variable_set, return all possible pair of global allocation?
//when calculate the permutation, can we store it and reuse to save time.






//we need to see if we can get the cost of each instruction directly from this function
//I hope it is not hard, but I am not sure.
//static float instruction_cost(i_assignment_ps &a);

static std::vector<var_t> unionVectors(std::vector<var_t>& vec1, 
                         std::vector<var_t>& vec2) 
{ 
    std::vector<var_t> ans; 
    // Declare the set to store the unique elements 
    std::set<var_t> s; 
    // insert elements from vector 1 into the set 
    for (int i = 0; i < vec1.size(); i++) { 
        s.insert(vec1[i]); 
    } 
    // insert elements from vector 2 into the set 
    for (int i = 0; i < vec2.size(); i++) { 
        s.insert(vec2[i]); 
    } 
    // Store the union of both the vectors into a resultant 
    // vector 
    for (auto it = s.begin(); it != s.end(); it++) { 
        ans.emplace_back(*it); 
    } 
    std::sort(ans.begin(),ans.end());
    return ans; 
} 


static void initial_after(cfg_t &cfg){
   boost::graph_traits<cfg_t>::vertex_iterator vi, vi_end;
   for (boost::tie(vi, vi_end) = vertices(cfg); vi != vi_end; ++vi){
          boost::graph_traits<cfg_t>::out_edge_iterator eout, eout_end;
            for (boost::tie(eout, eout_end) = out_edges(*vi, cfg); eout != eout_end; ++eout){
               vertex target=boost::target(*eout, cfg);
               cfg[*vi].after=unionVectors(cfg[*vi].after,cfg[target].alive);
            }
            std::sort(cfg[*vi].after.begin(),cfg[*vi].after.end());

   }

}

static void calcSubset(f& A, std::vector<f>& res,
                f& subset, int index)
{
    // Add the current subset to the result list
    res.push_back(subset);
 
    // Generate subsets by recursively including and
    // excluding elements
    for (int i = index; i < A.size(); i++) {
        // Include the current element in the subset
        subset.push_back(A[i]);
 
        // Recursively generate subsets with the current
        // element included
        calcSubset(A, res, subset, i + 1);
 
        // Exclude the current element from the subset
        // (backtracking)
        subset.pop_back();
    }
}

static std::vector<f> generate_permutation(f variables){
   std::vector<f> results;
   results.push_back(variables);
   while (std::next_permutation(variables.begin(), variables.end())){

      results.push_back(variables);}
   
   
   return results;
}
 
static bool compare(f a, f b){
   if (a.size()>b.size()){
      return false;
   }else if(a.size()<b.size()){
      return true;
   }else{
      for(int i=0;i<a.size();++i){
         if(a[i]>b[i]){
            return false;
         }else if(a[i]<b[i]){
            return true;
         }
      }
   
   }
}

std::vector<f > subsets(f& A)
{
   f subset;
    std::vector<f > res;
    int index = 0;
    calcSubset(A, res, subset, index);
    std::sort(res.begin(),res.end(),compare);
    return res;
}

static f extend_glob(f v_origin, f v_new, f origin_glob){
   for(auto i:v_new){
    if (origin_glob[i]==-3)
      {
         origin_glob[i]=-1;
      }
   }
   return origin_glob;
}

static void generate_possibility(f variables,int n){
  // std::cout<<"begin generate_possibility"<<std::endl;
 //  std::vector<f> results;
   std::vector<f> sub_set=subsets(variables);
   for(auto sub:sub_set){
      std::vector<f> globs;
      if(sub.size()<=MAX_NUM_REGS){
        f v;
        for(int i=0; i<MAX_NUM_REGS;++i){
          v.emplace_back(-1);
        }
        int len=sub.size();
        for(int i=0;i<len;++i){
          v[MAX_NUM_REGS-len+i]=sub[i];
        }
        std::vector<f> p=generate_permutation(v);
         for (auto i:p){
         f global;
         convert_to_global(i,sub,global,n);
         globs.emplace_back(global);
        }
      }
        std::vector<f> sub_sub_set=subsets(sub);
        for(auto i:sub_sub_set){
         if (i!=sub){
            for(auto j:permutation_map[i]){
               f g=extend_glob(i,sub,j);
               globs.emplace_back(g);
            }
         }
        }
        permutation_map[sub]=globs;
        
      //  results.reserve(results.size() + distance(p.begin(),p.end()));
       // results.insert(results.end(),p.begin(),p.end());
        //results.push_back(v);
      
   }
 //  std::cout<<"finish generate_possibility"<<std::endl;
}

static f get_partial_global(f global, f variables){
 //  std::cout<<"begin get_partial_global"<<std::endl;
   f result;
   result.resize(global.size(),-3);
   for(auto i:variables){
      result[i]=global[i];
   }
  // std::cout<<"finish get_partial_global"<<std::endl;
   return result;
}

//this function is used to combine two assignment_ps_list while series merge
static void combine_assignment_ps_list_series(ps_cfg_t &a, ps_cfg_t &b, ps_cfg_t &c){
   //assignment_ps first_a=a.assignments.begin()->second;
   //assignment_ps first_b=b.assignments.begin()->second;
   assignment_ps_map::iterator ita=a.assignments.begin();
   assignment_ps_map::iterator ita_end=a.assignments.end();
   if (a.variables==b.variables){
      for (;ita!=ita_end;++ita){
         c.assignments.emplace(std::make_pair(ita->first,assignment_ps(ita->second.s+b.assignments[ita->first].s,ita->second.begin_cost,b.assignments[ita->first].end_cost)));
      }
         return;

   }
   //std::cout<<"series size1 : "<<c.size()<<std::endl;
   f v_n=c.variables;
   assignment_ps_map::iterator itb;
   assignment_ps_map::iterator itb_begin=b.assignments.begin();
   assignment_ps_map::iterator itb_end=b.assignments.end();
   assignment_ps_map::iterator itc;
   for (;ita!=ita_end;++ita){
      for(itb=itb_begin;itb!=itb_end;++itb){
         f new_g;
         if_f_match(ita->first,itb->first,new_g);
         if(new_g[0]==-2){
            continue;
         }
         new_g=get_partial_global(new_g,v_n);
         float s=ita->second.s+itb->second.s;
         itc=c.assignments.find(new_g);
         if(itc==c.assignments.end()||itc->second.s>s){
         //if(c.find(new_g)==c.end()||c[new_g].s>s){
            c.assignments[new_g]=assignment_ps(s,ita->second.begin_cost,itb->second.end_cost);
         }
      }
   }
    //  std::cout<<"series size2 : "<<c.size()<<std::endl;
}

static void combine_assignment_ps_list_parallel(ps_cfg_t &a, ps_cfg_t &b, ps_cfg_t &c){
//assignment_ps first_a=a.begin()->second;
   assignment_ps_map::iterator ita=a.assignments.begin();
   assignment_ps_map::iterator ita_end=a.assignments.end();
   for (;ita!=ita_end;++ita){
     //float s=a[i].s+b[i].s-a[i].end_cost-a[i].begin_cost;
       c.assignments.emplace(std::make_pair(ita->first,assignment_ps(ita->second.s+b.assignments[ita->first].s-ita->second.begin_cost-ita->second.end_cost,ita->second.begin_cost,ita->second.end_cost)));
     
   }   
}

static void combine_assignment_ps_list_loop(ps_cfg_t &a, ps_cfg_t &b, ps_cfg_t &c){
  // std::cout<<"begin combine_assignment_ps_list_loop"<<std::endl;
   assignment_ps_map::iterator ita=b.assignments.begin();
   assignment_ps_map::iterator ita_end=b.assignments.end();

   if (b.variables==a.variables){
   for (;ita!=ita_end;++ita){
       c.assignments.emplace(std::make_pair(ita->first,assignment_ps(ita->second.s+a.assignments[ita->first].s,ita->second.begin_cost,ita->second.end_cost)));
   }
  // std::cout<<"loop size1 : "<<c.size()<<std::endl;
   return;
}
 f v_n=c.variables;
assignment_ps_map::iterator itb;
assignment_ps_map::iterator itb_begin=a.assignments.begin();
   assignment_ps_map::iterator itb_end=a.assignments.end();
assignment_ps_map::iterator itc;
for (;ita!=ita_end;++ita){
   for(itb=itb_begin;itb!=itb_end;++itb){
      f new_g;
      if_f_match(ita->first,itb->first,new_g);
      if(new_g[0]==-2){
         continue;
      }
      float s=ita->second.s+itb->second.s;
      new_g=get_partial_global(new_g,v_n);
      itc=c.assignments.find(new_g);
      if(itc==c.assignments.end()||itc->second.s>s){
      //if(c.find(new_g)==c.end()||c[new_g].s>s){
         c.assignments[new_g]=assignment_ps(s,ita->second.begin_cost,ita->second.end_cost);
      }
   }
   }
  // std::cout<<"loop size2 : "<<c.size()<<std::endl;
  // std::cout<<"finish combine_assignment_ps_list_loop"<<std::endl;
}

template <class I_t>
static float instruction_cost_easy(const f & global, cfg_node &node, const I_t &I);



static std::vector<f> generate_p_w(f variables){
   std::vector<f> results;
   f v;
   f current;
   v.resize(MAX_NUM_REGS,-1);
   results.push_back(v);
   for(auto i:variables){
      current.push_back(i);
      std::vector<f> new_results;
      for(auto j:results){
         for(int k=0;k<MAX_NUM_REGS;++k){
            if(j[k]==-1){
               f newf=j;
               newf[k]=i;
               new_results.push_back(newf);
            }

      }
   }
   results.insert(results.end(),new_results.begin(),new_results.end());
   }
 return results;     
}


template <class I_t>
static void initlize_assignment_ps_list(ps_cfg_t &a, I_t &I){
    
   std::vector<f>  begin_p=generate_p_w(a.begin_v);

   cfg_node node=cfg_map[a.index][0];
   int n = boost::num_vertices(I);
   for(auto i:begin_p){
         f global;
         convert_to_global(i,a.begin_v,global,n);

         a.assignments.emplace(std::make_pair(global, assignment_ps(instruction_cost_easy(global,node,I))));
   }
  // std::cout<<"basic block size: "<< a.assignments.size()<<std::endl;
}


template <class I_t>
static void generate_spcfg(ps_cfg_t &ps_cfg, I_t &I){
 // if (ps_cfg.make_series){
  //    series_addition(ps_cfg);
   //   return;
   //}

   if (ps_cfg.assignments.size() == 0){
       if(ps_cfg.left==-1 || ps_cfg.right==-1){
        // std::cout<<"1"<<std::endl;
         initlize_assignment_ps_list(ps_cfg, I);
         return;
      }

      if (ps_cfg_map[ps_cfg.left].assignments.size() == 0){
       //  std::cout<<"2"<<std::endl;
         generate_spcfg(ps_cfg_map[ps_cfg.left],I);
      }
      if (ps_cfg_map[ps_cfg.right].assignments.size() == 0){
       //  std::cout<<"3"<<std::endl;
         generate_spcfg(ps_cfg_map[ps_cfg.right],I);
      }
      switch (ps_cfg.type){
         case 1:
         //  std::cout<<"4"<<std::endl;
            combine_assignment_ps_list_series(ps_cfg_map[ps_cfg.left], ps_cfg_map[ps_cfg.right],ps_cfg);
         //  std::cout<<"current optimal:"<<get_optimal(ps_cfg,I).s<<std::endl;
            break;
         case 2:
         //   std::cout<<"5"<<std::endl;
            combine_assignment_ps_list_parallel(ps_cfg_map[ps_cfg.left], ps_cfg_map[ps_cfg.right],ps_cfg);
        //   std::cout<<"current optimal:"<<get_optimal(ps_cfg,I).s<<std::endl;
            break;
         case 3:
         //   std::cout<<"6"<<std::endl;
            combine_assignment_ps_list_loop(ps_cfg_map[ps_cfg.left], ps_cfg_map[ps_cfg.right],ps_cfg);
         //  std::cout<<"current optimal:"<<get_optimal(ps_cfg,I).s<<std::endl;
            break;
         default:
            break;
      }
   }
}

static assignment_ps get_optimal(ps_cfg_t &ps_cfg){
   assignment_ps_map a = ps_cfg.assignments;
   assignment_ps b;
   b.s = std::numeric_limits<float>::infinity();
   f c;
   for(auto i:a){
      //std::cout<<"i.second.s:"<<i.second.s<<std::endl;
      //std::cout << std::endl;
      if(b.s > i.second.s){
         b = i.second;
         c=i.first;
      }
   }

   return b;
}