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
struct i_assignment_ps{
   f registers_begin;
   f global_regs;
   cfg_node *node; //the corresponding node(with ic) in the cfg
   float cost; //cost of the assignment

   i_assignment_ps(){
   // std::cout<<"i_assignment_ps constructor"<<std::endl;
      registers_begin.clear();
      for(int i=0; i<MAX_NUM_REGS; i++){
         registers_begin.push_back(-1);
      }
     // std::cout<<"i_assignment_ps constructor end"<<std::endl;
      node = NULL;
      cost = std::numeric_limits<float>::infinity();
   }

   
};

//assignment for graph
struct assignment_ps{
   float s; //cost
   //std::vector<i_assignment_ps> insts; //assignments for each instruction
   float begin_cost;
    float end_cost;
  // f global_regs;

   assignment_ps(){
    //std::cout<<"assignment_ps constructor"<<std::endl;
      s = std::numeric_limits<float>::infinity();
      //insts.clear();
     // std::cout<<"assignment_ps constructor end"<<std::endl;
   }
};


typedef std::map<std::pair<f,f>,assignment_ps> assignment_ps_map;



struct ps_cfg_t{

  vertex begin; //s node
  vertex end; //t node
  cfg_t* cfg; //graph contains s and t node
  merge_type type;
  int left;
  int right;
  int index;
  int parent;
  varset_t begin_v;
  varset_t end_v;
  assignment_ps_map assignments;
};

std::map<int,ps_cfg_t> ps_cfg_map;
int ps_cfg_count=0;
std::map<int,cfg_t> cfg_map;
int cfg_count=0;

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



static ps_cfg_t init_ps_cfg(cfg_t &cfg, vertex begin_node, vertex end_node, int index, int parent){
  ps_cfg_t ps_cfg;
  ps_cfg.cfg=&cfg;
  ps_cfg.begin=begin_node;
  ps_cfg.end=end_node;
  ps_cfg.type=0;
  ps_cfg.left=-1;
  ps_cfg.right=-1;
  ps_cfg.index=index;
  ps_cfg.parent=parent;
  ps_cfg.assignments.clear();
  ps_cfg.begin_v=cfg[begin_node].alive;
  ps_cfg.end_v=cfg[end_node].after;
  
    return ps_cfg;
}

//break the graph into two series part
static void break_graph_series(cfg_t &cfg, vertex begin_node, vertex end_node, cfg_t &cfg_1, cfg_t &cfg_2,ps_cfg_t &ps_cfg){
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
  cfg_map[cfg_count]=cfg_1;
  ps_cfg_t right=init_ps_cfg(cfg_map[cfg_count],find_vertex_from_node(cfg[end_node],cfg_1), find_vertex_from_node(cfg[*(vi_end-1)],cfg_1),ps_cfg_count,ps_cfg.index);
  cfg_count++;
  ps_cfg_map[ps_cfg_count]=right;
  ps_cfg.right=ps_cfg_count;
  ps_cfg_count++;
  f=vi;
  vi=next;
  for (next = vi; next != vi_end; vi=next) {
      ++next;
      boost::clear_vertex(find_vertex_from_node(cfg[*vi], cfg_2),cfg_2);
      boost::remove_vertex(find_vertex_from_node(cfg[*vi], cfg_2),cfg_2);
  }
  cfg_map[cfg_count]=cfg_2;
  ps_cfg_t left=init_ps_cfg(cfg_map[cfg_count],find_vertex_from_node(cfg[begin_node], cfg_2), find_vertex_from_node(cfg[*f],cfg_2),ps_cfg_count,ps_cfg.index);
  cfg_count++;
  ps_cfg_map[ps_cfg_count]=left;
  ps_cfg.left=ps_cfg_count;
  ps_cfg_count++;
}

static void break_graph_parallel(cfg_t &cfg, vertex begin_node, vertex end_node, cfg_t &cfg_1, cfg_t &cfg_2,ps_cfg_t &ps_cfg){
  boost::graph_traits<cfg_t>::vertex_iterator vi, vi_end, next;
  boost::tie(vi, vi_end) = boost::vertices(cfg);
  ps_cfg.type=2;
  boost::graph_traits < cfg_t >::out_edge_iterator ei, ei_end;
  vertex p_s=0;
  for (boost::tie(ei, ei_end) = boost::out_edges(begin_node, cfg); ei != ei_end; ++ei) {
    vertex target = boost::target ( *ei, cfg );
    if (p_s==0 || target>p_s){
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
  cfg_map[cfg_count]=cfg_1;
  ps_cfg_t left=init_ps_cfg(cfg_map[cfg_count],find_vertex_from_node(cfg[begin_node], cfg_1), find_vertex_from_node(cfg[end_node], cfg_1),ps_cfg_count,ps_cfg.index);
  cfg_count++;
  ps_cfg_map[ps_cfg_count]=left;
  ps_cfg.left=ps_cfg_count;
  ps_cfg_count++;
  vi=next;
  --vi_end;
  for (next = vi; next != vi_end; vi=next) {
      ++next;
      boost::clear_vertex(find_vertex_from_node(cfg[*vi], cfg_2),cfg_2);
      boost::remove_vertex(find_vertex_from_node(cfg[*vi], cfg_2),cfg_2);
  }
  cfg_map[cfg_count]=cfg_2;
  ps_cfg_t right=init_ps_cfg(cfg_map[cfg_count],find_vertex_from_node(cfg[begin_node], cfg_2), find_vertex_from_node(cfg[end_node], cfg_2),ps_cfg_count,ps_cfg.index);
  cfg_count++;
  ps_cfg_map[ps_cfg_count]=right;
  ps_cfg.right=ps_cfg_count;
  ps_cfg_count++;
}

static vertex find_parallel_end(cfg_t &cfg){
  boost::graph_traits<cfg_t>::vertex_iterator vi, vi_end, next;
  boost::tie(vi, vi_end) = boost::vertices(cfg);
  int count=0;
  for (next = vi; next != vi_end; vi=next) {
    ++next;
    if (boost::in_degree(*vi,cfg)==2 ){
         boost::graph_traits < cfg_t >::in_edge_iterator ei, ei_end;
         int flag=0;
         for (boost::tie(ei, ei_end) = boost::in_edges(*vi, cfg); ei != ei_end; ++ei) {
             vertex source = boost::source ( *ei, cfg );

             if(source < *vi){
               flag++;
             }
         }
         if (flag==2){
            count=count-1;
            if (count==0){
              return *vi;
            }
         }else{
            next=next+1;
            boost::graph_traits < cfg_t >::out_edge_iterator eei, eei_end;
            vertex l_end=0;
            for (boost::tie(eei, eei_end) = boost::out_edges(*next, cfg); eei != eei_end; ++eei) {
              vertex target = boost::target ( *eei, cfg );
              if( target>l_end){
                l_end=target;
              }
            }
            while (*next!=l_end){
              next++;
            }
          }
         
    }
    if(boost::out_degree(*vi,cfg)==2){
      count++;
    }
  }
  return -1;
}

static void break_graph_loop(cfg_t &cfg, vertex begin_node, vertex end_node, cfg_t &cfg_1, cfg_t &cfg_2,ps_cfg_t &ps_cfg){
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
  cfg_map[cfg_count]=cfg_1;
  ps_cfg_t left=init_ps_cfg(cfg_map[cfg_count],find_vertex_from_node(cfg[begin_node],cfg_1), find_vertex_from_node(cfg[*vi],cfg_1),ps_cfg_count,ps_cfg.index);
  cfg_count++;
  ps_cfg_map[ps_cfg_count]=left;
  ps_cfg.left=ps_cfg_count;
  ps_cfg_count++;
  vi=f;
  for(next=vi; next!=vi_end; vi=next){
    ++next;
    boost::clear_vertex(find_vertex_from_node(cfg[*vi], cfg_2),cfg_2);
    boost::remove_vertex(find_vertex_from_node(cfg[*vi], cfg_2),cfg_2);
  }
  cfg_map[cfg_count]=cfg_2;
  ps_cfg_t right=init_ps_cfg(cfg_map[cfg_count],find_vertex_from_node(cfg[*f],cfg_2), find_vertex_from_node(cfg[*(vi_end-1)],cfg_2),ps_cfg_count,ps_cfg.index);
  cfg_count++;
  ps_cfg_map[ps_cfg_count]=right;
  ps_cfg.right=ps_cfg_count;
  ps_cfg_count++;
}

static void convert_cfg_to_spcfg_one_step(ps_cfg_t &pscfg){
  //convert the original cfg to the root node of ps_cfg
  cfg_t cfg=*(pscfg.cfg);
  boost::graph_traits<cfg_t>::vertex_iterator vi, vi_end, next;
  boost::tie(vi, vi_end) = boost::vertices(cfg);
  //pscfg=init_ps_cfg(cfg, *vi, *(vi_end-1));
  if (*vi == *(vi_end-1)) {
    return ;
  }

  next=vi;
  next++;
  //series merge - boost::in_degree=1, boost::out_degree=1
  if (boost::in_degree(*vi,cfg)==0 && boost::out_degree(*vi, cfg) == 1){
   // std::cout<<"series merge break!"<<std::endl;
    //boost::write_graphviz(std::cout, cfg,boost::make_label_writer(boost::get(&cfg_node::ic,cfg )));
    cfg_t cfg_1, cfg_2;
    boost::copy_graph(cfg, cfg_1);
    boost::copy_graph(cfg, cfg_2);
    break_graph_series(cfg, *vi, *next, cfg_1, cfg_2, pscfg);

    //boost::write_graphviz(std::cout,cfg_1);
    return;
  }
  //start of parallel merge - in degree=1, out degree=2, both next verdexes with larger index
  //end of parallel merge - in degree=2, out degree=1, both previous vertexes with smaller index
  if (boost::in_degree(*vi,cfg)==0 && boost::out_degree(*vi,cfg)==2){
   // std::cout<<"parallel merge break!"<<std::endl;
    vertex p_end=find_parallel_end(cfg);
    if(*(vi_end-1)!=p_end){
      //if the graph is built by parallel merge then series merge
      cfg_t cfg_1, cfg_2;
      boost::copy_graph(cfg, cfg_1);
      boost::copy_graph(cfg, cfg_2);
      next=vi;
      while(*next!=p_end){
        next++;
      }
      break_graph_series(cfg, *vi, *(next+1), cfg_1, cfg_2, pscfg);
      return ;
    }else{
      //if the graph is built by parallel merge directly
      cfg_t cfg_1, cfg_2;
      boost::copy_graph(cfg, cfg_1);
      boost::copy_graph(cfg, cfg_2);
      break_graph_parallel(cfg, *vi, p_end, cfg_1, cfg_2, pscfg);
    }
    return ;
  }
  //start of loop merge - in degree=2, out degree=1, one previous vertex is larger and one is smaller
  //end of loop merge - in degree=1, out degree=2, both next vertexes with larger index.
  if (boost::in_degree(*vi,cfg)==1){
    if (boost::out_degree(*(vi+2),cfg)==1){
      //std::cout<<"loop merge break!"<<std::endl;
      cfg_t cfg_1, cfg_2;
      boost::copy_graph(cfg, cfg_1);
      boost::copy_graph(cfg, cfg_2);
      break_graph_loop(cfg, *vi, *(vi+2), cfg_1, cfg_2, pscfg);

      return ;
    }else{
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
      break_graph_series(cfg, *vi, l_end, cfg_1, cfg_2, pscfg);
      return ;
    }
  }
  return ;
}

static void convert_cfg_to_spcfg(ps_cfg_t &pscfg){
  //ps_cfg_map.clear();
  //cfg_map.clear();
  //ps_cfg_count=0;
  //cfg_count=0;
  convert_cfg_to_spcfg_one_step(pscfg);
  if(pscfg.left!=-1 || pscfg.right!=-1){
    if (pscfg.left!=-1){
      convert_cfg_to_spcfg(ps_cfg_map[pscfg.left]);
    }
    if (pscfg.right!=-1){
      convert_cfg_to_spcfg(ps_cfg_map[pscfg.right]);
    }
  }
}