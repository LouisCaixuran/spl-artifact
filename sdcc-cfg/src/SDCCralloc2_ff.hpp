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

static bool if_f_match(f f1,f f2){
   std::set<short int> v_not_in_1;
   std::set<short int> v_not_in_2;
   for(int i=0;i<MAX_NUM_REGS;i++){
      if(f1[i]!=f2[i]){
         if(f1[i]!=-1){
            if (v_not_in_1.find(f1[i]) != v_not_in_1.end()){
               return false;
            }
            v_not_in_2.insert(f1[i]);
         }
         if(f2[i]!=-1){
            if (v_not_in_2.find(f2[i]) != v_not_in_2.end()){
               return false;
            }
            v_not_in_1.insert(f2[i]);
         }
      }
   }
   return true;
}


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
        ans.push_back(*it); 
    } 
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
      results.push_back(variables);
   }
   
   return results;
}
 
std::vector<f > subsets(f& A)
{
   f subset;
    std::vector<f > res;
    int index = 0;
    calcSubset(A, res, subset, index);
    return res;
}

static std::vector<f> generate_possibility(f variables){
  // std::cout<<"begin generate_possibility"<<std::endl;
   std::vector<f> results;
   std::vector<f> sub_set=subsets(variables);
   for(auto sub:sub_set){
      if(sub.size()<=MAX_NUM_REGS){
        f v;
        for(int i=0; i<MAX_NUM_REGS;i++){
          v.push_back(-1);
        }
        int len=sub.size();
        for(int i=0;i<len;i++){
          v[MAX_NUM_REGS-len+i]=sub[i];
        }

        std::vector<f> p=generate_permutation(v);

        results.reserve(results.size() + distance(p.begin(),p.end()));
        results.insert(results.end(),p.begin(),p.end());
        //results.push_back(v);
      }
   }
 //  std::cout<<"finish generate_possibility"<<std::endl;

   return results;
}



//this function is used to combine two assignment_ps_list while series merge
static assignment_ps_map combine_assignment_ps_list_series(ps_cfg_t a, ps_cfg_t b){
   assignment_ps_map c;
  for(auto i :a.assignments){
   f beginS=i.first.first;
   f endS=i.first.second;
   for(auto j:b.assignments){
      f beginT=j.first.first;
      f endT=j.first.second;
      if (!if_f_match(endS,beginT)){
         continue;
      }
      assignment_ps ac;
      ac.s=i.second.s+j.second.s;
      ac.begin_cost=i.second.begin_cost;
      ac.end_cost=j.second.end_cost;
       if(c[std::pair<f,f>(beginS,endT)].s > ac.s){
            //   ac.global_regs=aa.global_regs;
              c[std::pair<f,f>(beginS,endT)] = ac;
      }

   }
  }
   //std::cout<<"combine_assignment_ps_list_series.size"<<c.size() <<std::endl;
   return c;
}

static assignment_ps_map combine_assignment_ps_list_parallel(ps_cfg_t a, ps_cfg_t b){
   assignment_ps_map c;

 for(auto i :a.assignments){
   f beginS=i.first.first;
   f endS=i.first.second;
     assignment_ps ab=b.assignments[std::pair<f,f>(beginS,endS)];
      if(ab.s==std::numeric_limits<float>::infinity()){
         continue;
      }
   assignment_ps ac;
   ac.s=i.second.s+ab.s-ab.end_cost-ab.begin_cost;
   ac.begin_cost=i.second.begin_cost;
   ac.end_cost=i.second.end_cost;
   c[std::pair<f,f>(beginS,endS)] = ac;
  }
   return c;
   
}

static assignment_ps_map combine_assignment_ps_list_loop(ps_cfg_t a, ps_cfg_t b){
  // std::cout<<"begin combine_assignment_ps_list_loop"<<std::endl;
   assignment_ps_map c;
   for(auto i : b.assignments){
      f beginS=i.first.first;
      f endS=i.first.second;
      assignment_ps ab=a.assignments[std::pair<f,f>(endS,beginS)];
      if(ab.s==std::numeric_limits<float>::infinity()){
         continue;
      }
      assignment_ps ac;
      ac.s=i.second.s+ab.s;
      ac.begin_cost=i.second.begin_cost;
      ac.end_cost=i.second.end_cost;
      c[std::pair<f,f>(beginS,endS)] = ac;
   }
   return c;
}

template <class I_t>
static float instruction_cost_easy(const i_assignment_ps &ia, cfg_node &node, const I_t &I);

static void initlize_assignment_ps_list(ps_cfg_t &a){
   assignment_ps_map c;

   std::vector<f> begin=generate_possibility(a.begin_v);
   //std::cout<<"begin size:"<<begin.size()<<std::endl;
   //std::cout<<"end size:"<<end.size()<<std::endl;
   //std::cout<<"finish generating"<<std::endl;
   for(auto i:begin){
         //std::cout<<"begin to get cost"<<std::endl;
         assignment_ps aa=assignment_ps();
         //std::cout<<"finish initial assignment_ps"<<std::endl;
         i_assignment_ps as=i_assignment_ps();
         // std::cout<<"finish initial assignment"<<std::endl;
         as.registers_begin = i;
         //std::cout<<"try to get node"<<std::endl;
         as.node=&((*(a.cfg))[a.begin]);
         as.global_regs.reserve(a.begin_v.size());
         for(auto i : a.begin_v){
            as.global_regs[i]=getIndex(as.registers_begin,i);
         }
         //std::cout<<"try to get cost"<<std::endl;
         as.cost = instruction_cost_easy(as,*(as.node),I);
         aa.s = as.cost;
         aa.begin_cost=as.cost;
         aa.end_cost=as.cost;
         //aa.begin_i = as;
         //aa.end_i = as;
         c[std::pair<f,f>(i,i)] = aa;
   }
  // std::cout<<"c.size()"<<c.size()<<std::endl;

   a.assignments = c;
}

static void generate_spcfg(ps_cfg_t &ps_cfg){

    if(ps_cfg.left==-1 || ps_cfg.right==-1){
        // std::cout<<"1"<<std::endl;
         initlize_assignment_ps_list(ps_cfg, I);
         return;
      }

   if (ps_cfg.assignments.size() == 0){
      if (ps_cfg_map[ps_cfg.left].assignments.size() == 0){
       //  std::cout<<"2"<<std::endl;
         generate_spcfg(ps_cfg_map[ps_cfg.left]);
      }
      if (ps_cfg_map[ps_cfg.right].assignments.size() == 0){
       //  std::cout<<"3"<<std::endl;
         generate_spcfg(ps_cfg_map[ps_cfg.right]);
      }
      switch (ps_cfg.type){
         case 1:
         //  std::cout<<"4"<<std::endl;
            ps_cfg.assignments = combine_assignment_ps_list_series(ps_cfg_map[ps_cfg.left], ps_cfg_map[ps_cfg.right]);
         //  std::cout<<"current optimal:"<<get_optimal(ps_cfg,I).s<<std::endl;
            break;
         case 2:
         //   std::cout<<"5"<<std::endl;
            ps_cfg.assignments = combine_assignment_ps_list_parallel(ps_cfg_map[ps_cfg.left], ps_cfg_map[ps_cfg.right]);
        //   std::cout<<"current optimal:"<<get_optimal(ps_cfg,I).s<<std::endl;
            break;
         case 3:
         //   std::cout<<"6"<<std::endl;
            ps_cfg.assignments = combine_assignment_ps_list_loop(ps_cfg_map[ps_cfg.left], ps_cfg_map[ps_cfg.right]);
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
   for(auto i:a){
      //std::cout<<"i.second.s:"<<i.second.s<<std::endl;
      //std::cout << std::endl;
      if(b.s > i.second.s){
         b = i.second;
      }
   }
   return b;
}