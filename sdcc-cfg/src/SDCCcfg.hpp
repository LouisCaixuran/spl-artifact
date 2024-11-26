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

extern "C"
{
#include "SDCCbtree.h"
}


typedef int merge_type;
typedef boost::graph_traits<cfg_t>::vertex_descriptor vertex;
typedef std::vector<var_t> f;

#define MAX_NUM_REGS 3

//assignment for one instuction


//assignment for graph
struct assignment_ps{
   float s; //cost
   //std::vector<i_assignment_ps> insts; //assignments for each instruction
   float begin_cost;
   float end_cost;
   
   assignment_ps(){
    //std::cout<<"assignment_ps constructor"<<std::endl;
      s = std::numeric_limits<float>::infinity();
      //insts.clear();
     // std::cout<<"assignment_ps constructor end"<<std::endl;
   }

   assignment_ps(float cost){
      s = cost;
      begin_cost = cost;
      end_cost = cost;
 

   }

   assignment_ps(float cost, float begin, float end){
      s = cost;
      begin_cost = begin;
      end_cost = end;

   }
};

typedef std::map<f,assignment_ps> assignment_ps_map;

std::map<f,std::vector<f>> permutation_map;


struct ps_cfg_t{

  vertex begin; //s node
  vertex end; //t node
  //cfg_t* cfg; //graph contains s and t node
  merge_type type;
  int left;
  int right;
  int index;
  int parent;
  varset_t begin_v;
  varset_t end_v;
  varset_t variables;
  assignment_ps_map assignments;

};

std::vector<ps_cfg_t> ps_cfg_map;
std::vector<cfg_t> cfg_map;
int cfg_count=1;

  static ps_cfg_t init_ps_cfg( vertex begin_node, vertex end_node, int index, int parent){
  ps_cfg_t pscfg;
  pscfg.begin=begin_node;
  pscfg.end=end_node;
  pscfg.type=0;
  pscfg.index=index;
  pscfg.parent=parent;
  pscfg.begin_v=cfg_map[index][begin_node].alive;
  pscfg.end_v=cfg_map[index][end_node].alive;
  std::set_union(pscfg.begin_v.begin(),pscfg.begin_v.end(),pscfg.end_v.begin(),pscfg.end_v.end(),std::inserter(pscfg.variables,pscfg.variables.begin()));
return pscfg;
}

static vertex find_vertex_from_node(cfg_node node,cfg_t &cfg){
  boost::graph_traits<cfg_t>::vertex_iterator vi, vi_end;
  boost::tie(vi, vi_end) = boost::vertices(cfg);
  for (; vi != vi_end; ++vi) {
    if (cfg[*vi].ic == node.ic) {
      return *vi;
    }
  }
  return 0;
}


//break the graph into two series part
static void break_graph_series( vertex begin_node, vertex end_node, cfg_t &cfg_1, cfg_t &cfg_2,ps_cfg_t &ps_cfg){
  //std::cout<<"break_graph_series begin"<<std::endl;
  cfg_t cfg=cfg_map[ps_cfg.index];
  boost::graph_traits<cfg_t>::vertex_iterator vi, vi_end, next,f;
  boost::tie(vi, vi_end) = boost::vertices(cfg);
  ps_cfg.type=1;

  for (next = vi; next != vi_end; vi=next){
      ++next;
      boost::clear_vertex(find_vertex_from_node(cfg[*vi], cfg_1),cfg_1);
      boost::remove_vertex(find_vertex_from_node(cfg[*vi], cfg_1),cfg_1);
      if (*next==end_node){
        break;
      }
  }
  cfg_map.push_back(cfg_1);
  ps_cfg_t right=init_ps_cfg(find_vertex_from_node(cfg[end_node],cfg_1), find_vertex_from_node(cfg[*(vi_end-1)],cfg_1),cfg_count,ps_cfg.index);
  ps_cfg.right=cfg_count;
  cfg_count++;
  ps_cfg_map.push_back(right);
  f=vi;
  vi=next;
  for (next = vi; next != vi_end; vi=next) {
      ++next;
      boost::clear_vertex(find_vertex_from_node(cfg[*vi], cfg_2),cfg_2);
      boost::remove_vertex(find_vertex_from_node(cfg[*vi], cfg_2),cfg_2);
  }
  cfg_map.push_back( cfg_2);
  ps_cfg_t left=init_ps_cfg(find_vertex_from_node(cfg[begin_node],cfg_2), find_vertex_from_node(cfg[*f],cfg_2),cfg_count,ps_cfg.index);
  ps_cfg.left=cfg_count;
  cfg_count++;
  ps_cfg_map.push_back(left);
//  std::cout<<"break_graph_series end"<<std::endl;
}

static void break_graph_parallel( vertex begin_node, vertex end_node, cfg_t &cfg_1, cfg_t &cfg_2,ps_cfg_t &ps_cfg){
 // std::cout<<"break_graph_parallel begin"<<std::endl;
 cfg_t cfg=cfg_map[ps_cfg.index];
  boost::graph_traits<cfg_t>::vertex_iterator vi, vi_end, next;
  boost::tie(vi, vi_end) = boost::vertices(cfg);
  ps_cfg.type=2;
  boost::graph_traits < cfg_t >::out_edge_iterator ei, ei_end;
  vertex p_s=0;
  for (boost::tie(ei, ei_end) = boost::out_edges(begin_node, cfg); ei != ei_end; ++ei) {
    vertex target = boost::target ( *ei, cfg );
    if (p_s==begin_node || target>p_s){
      p_s=target;
    }
  }
  ++vi;
  for (next = vi; next != vi_end; vi=next) {
      ++next;
      boost::clear_vertex(find_vertex_from_node(cfg[*vi], cfg_1),cfg_1);
      boost::remove_vertex(find_vertex_from_node(cfg[*vi], cfg_1),cfg_1);
      if (*next==p_s){
        break;
      }
  }
  cfg_map.push_back(cfg_1);
  ps_cfg_t left=init_ps_cfg(find_vertex_from_node(cfg[begin_node],cfg_1), find_vertex_from_node(cfg[end_node],cfg_1) ,cfg_count,ps_cfg.index);
  ps_cfg.left=cfg_count;
  cfg_count++;
  ps_cfg_map.push_back(left);
  --vi_end;
  if(next==vi_end){
    boost::remove_edge(*vi,*vi_end,cfg_2);
  }else{
  vi=next;
  for (next = vi; next != vi_end; vi=next) {
      ++next;
      boost::clear_vertex(find_vertex_from_node(cfg[*vi], cfg_2),cfg_2);
      boost::remove_vertex(find_vertex_from_node(cfg[*vi], cfg_2),cfg_2);
  }}
  cfg_map.push_back(cfg_2);
  ps_cfg_t right=init_ps_cfg(find_vertex_from_node(cfg[begin_node],cfg_2), find_vertex_from_node(cfg[end_node],cfg_2) ,cfg_count,ps_cfg.index);
  ps_cfg.right=cfg_count;
  cfg_count++;
  ps_cfg_map.push_back(right);
//  std::cout<<"break_graph_parallel end"<<std::endl;

}

static vertex find_parallel_end(cfg_t &cfg){
  boost::graph_traits<cfg_t>::vertex_iterator vi, vi_end, next;
  boost::tie(vi, vi_end) = boost::vertices(cfg);
  vertex vertex_out=0;
  boost::graph_traits < cfg_t >::out_edge_iterator ei, ei_end;
  boost::tie(ei, ei_end) = boost::out_edges(*vi, cfg);
  for(;ei!=ei_end;++ei){
    vertex target = boost::target ( *ei, cfg );
    if (target>vertex_out){
      vertex_out=target;
    }
  }
  if(vertex_out==0){
    throw std::invalid_argument("invalid parallel");
  }
  vertex_out--;
  if(boost::out_degree(vertex_out,cfg)!=1){
    throw std::invalid_argument("invalid parallel");
  }
  boost::tie(ei, ei_end) = boost::out_edges(vertex_out, cfg);
  vertex target = boost::target ( *ei, cfg );
  return target;
  
}

static void break_graph_loop(vertex begin_node, vertex end_node, cfg_t &cfg_1, cfg_t &cfg_2,ps_cfg_t &ps_cfg){
//  std::cout<<"break_graph_loop begin"<<std::endl;
  cfg_t cfg=cfg_map[ps_cfg.index];
  ps_cfg.type=3;
  boost::graph_traits<cfg_t>::vertex_iterator vi, vi_end, next,f;
  boost::tie(vi, vi_end) = boost::vertices(cfg);
  f=vi+3;
  for(next=vi; next!=vi_end; vi=next){
    ++next;
    boost::clear_vertex(find_vertex_from_node(cfg[*vi], cfg_1),cfg_1);
    boost::remove_vertex(find_vertex_from_node(cfg[*vi], cfg_1),cfg_1);
    if (next==f){
      break;
    }
  }
  cfg_map.push_back(cfg_1);
  ps_cfg_t left=init_ps_cfg(find_vertex_from_node(cfg[*f],cfg_1), find_vertex_from_node(cfg[*(vi_end-1)],cfg_1),cfg_count,ps_cfg.index);
    ps_cfg.left=cfg_count;

  cfg_count++;
  ps_cfg_map.push_back(left);
  vi=f;
  for(next=vi; next!=vi_end; vi=next){
    ++next;
    boost::clear_vertex(find_vertex_from_node(cfg[*vi], cfg_2),cfg_2);
    boost::remove_vertex(find_vertex_from_node(cfg[*vi], cfg_2),cfg_2);
  }
  cfg_map.push_back(cfg_2);
  ps_cfg_t right=init_ps_cfg(find_vertex_from_node(cfg[begin_node],cfg_2), find_vertex_from_node(cfg[begin_node+2],cfg_2),cfg_count,ps_cfg.index);
  ps_cfg.right=cfg_count;

  cfg_count++;
  ps_cfg_map.push_back(right);
 // std::cout<<"break_graph_loop end"<<std::endl;
}

static void check_cfg(cfg_t cfg){
  boost::graph_traits<cfg_t>::vertex_iterator vi, vi_end;
  boost::tie(vi, vi_end) = boost::vertices(cfg);
  int num_begin=0;
  int num_end=0;
  int num_loop_or_parallel=0;
  for (; vi != vi_end; ++vi) {
    if(boost::in_degree(*vi,cfg)>2 || boost::out_degree(*vi,cfg)>2){
      throw std::invalid_argument("break or continue");
    }
   if (boost::in_degree(*vi,cfg)==0){
     num_begin++;
   }
    if (boost::out_degree(*vi,cfg)==0){
      num_end++;
    }
    if (boost::in_degree(*vi,cfg)==2){
      if( boost::out_degree(*vi,cfg)==1){
         num_loop_or_parallel++;
       //  std::cout<<"current node: n++ "<<*vi<<std::endl;
         }
      else{
        throw std::invalid_argument("invalid begin loop/end parallel");
      }
    }
   if (boost::out_degree(*vi,cfg)==2){
      if( boost::in_degree(*vi,cfg)==1){
         num_loop_or_parallel--;
    //              std::cout<<"current node: "<<*vi<<std::endl;
}
      else{
        throw std::invalid_argument("invalid end loop/begin parallel");
      }
    }
  }
  if(num_begin!=1 || num_end!=1 || num_loop_or_parallel!=0){
   // std::cout<<"num_begin: "<<num_begin<<std::endl;
  //  std::cout<<"num_end: "<<num_end<<std::endl;
  //  std::cout<<"num_loop_or_parallel: "<<num_loop_or_parallel<<std::endl;
    throw std::invalid_argument("invalid cfg");
  }

}

static void convert_cfg_to_spcfg_one_step(ps_cfg_t &pscfg){
  //convert the original cfg to the root node of ps_cfg
  cfg_t cfg=cfg_map[pscfg.index];
  //pscfg=init_ps_cfg(cfg, *vi, *(vi_end-1));
 // std::cout<<"get cfg with size: "<<boost::num_vertices(cfg)<<std::endl;
  if (boost::num_vertices(cfg)==1){ 
    //std::cout<<"basic block break!"<<std::endl;
    pscfg.left=-1;
    pscfg.right=-1;
  //  std::cout<<"basic block"<<std::endl;
    return ;
  }
  boost::graph_traits<cfg_t>::vertex_iterator vi, vi_end, next;
  boost::tie(vi, vi_end) = boost::vertices(cfg);

  next=vi;
  next++;
  //series merge - boost::in_degree=1, boost::out_degree=1
  if (boost::in_degree(*vi,cfg)==0 && boost::out_degree(*vi, cfg) == 1){
  //  std::cout<<"series merge break!"<<std::endl;
    //boost::write_graphviz(std::cout, cfg,boost::make_label_writer(boost::get(&cfg_node::ic,cfg )));
    cfg_t cfg_1, cfg_2;
    boost::copy_graph(cfg, cfg_1);
    boost::copy_graph(cfg, cfg_2);
    break_graph_series( *vi, *next, cfg_1, cfg_2, pscfg);

    //boost::write_graphviz(std::cout,cfg_1);
    return;
  }
  //start of parallel merge - in degree=1, out degree=2, both next verdexes with larger index
  //end of parallel merge - in degree=2, out degree=1, both previous vertexes with smaller index
  if (boost::in_degree(*vi,cfg)==0 && boost::out_degree(*vi,cfg)==2){
    vertex p_end=find_parallel_end(cfg);
    if(*(vi_end-1)!=p_end){
    //  std::cout<<"preparallel merge break!"<<std::endl;

      //if the graph is built by parallel merge then series merge
      cfg_t cfg_1, cfg_2;
      boost::copy_graph(cfg, cfg_1);
      boost::copy_graph(cfg, cfg_2);

      break_graph_series( *vi, p_end+1, cfg_1, cfg_2, pscfg);
      return ;
    }else{
      //if the graph is built by parallel merge directly
   //   std::cout<<"parallel merge break!"<<std::endl;
      cfg_t cfg_1, cfg_2;
      boost::copy_graph(cfg, cfg_1);
      boost::copy_graph(cfg, cfg_2);
      break_graph_parallel( *vi, p_end, cfg_1, cfg_2, pscfg);
    }
    return ;
  }
  //start of loop merge - in degree=2, out degree=1, one previous vertex is larger and one is smaller
  //end of loop merge - in degree=1, out degree=2, both next vertexes with larger index.
  if (boost::in_degree(*vi,cfg)==1){
    if (boost::out_degree(*(vi+2),cfg)==1){
  //    std::cout<<"loop merge break!"<<std::endl;
      cfg_t cfg_1, cfg_2;
      boost::copy_graph(cfg, cfg_1);
      boost::copy_graph(cfg, cfg_2);
      break_graph_loop( *vi, *(vi+2), cfg_1, cfg_2, pscfg);

      return ;
    }else{
      if((boost::out_degree(*(vi+2),cfg)!=2)){
        throw std::invalid_argument("invalid cfg while loop graph");
      }
  //    std::cout<<"preloop merge break!"<<std::endl;
      boost::graph_traits < cfg_t >::out_edge_iterator eei, eei_end;
      //std::cout<<"series merge break!2"<<std::endl;
      //boost::write_graphviz(std::cout, cfg,boost::make_label_writer(boost::get(&cfg_node::ic,cfg )));
      vertex l_end=0;
      next=vi+2;
      for (boost::tie(eei, eei_end) = boost::out_edges(*next, cfg); eei != eei_end; ++eei) {
        vertex target = boost::target ( *eei, cfg );
        //std::cout<<"target: "<<target<<std::endl;
        if( target>l_end){
          l_end=target;
        }
      }
      //std::cout<<"l_end: "<<l_end<<std::endl;
      cfg_t cfg_1, cfg_2;
      boost::copy_graph(cfg, cfg_1);
      boost::copy_graph(cfg, cfg_2);
      break_graph_series( *vi, l_end, cfg_1, cfg_2, pscfg);
      return ;
    }
  }
  throw std::invalid_argument("invalid cfg while breaking graph");
}

static void convert_cfg_to_spcfg(ps_cfg_t &pscfg){
  //ps_cfg_map.clear();
  //cfg_map.clear();
  //ps_cfg_count=0;
  //cfg_count=0;

  convert_cfg_to_spcfg_one_step(pscfg);
  int left=pscfg.left;
  int right=pscfg.right;
  if(left!=-1 || right!=-1){
    if (left!=-1){
      convert_cfg_to_spcfg(ps_cfg_map[left]);
    }
    if (right!=-1){
      convert_cfg_to_spcfg(ps_cfg_map[right]);
    }
  }

  //mark_series_pscfg(pscfg);
}

