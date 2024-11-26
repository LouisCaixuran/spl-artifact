#define TD_SALLOC
#define CH_SALLOC

#include "SDCCralloc.hpp"
#include "parallel_hashmap/phmap.h"
#include <boost/container/small_vector.hpp>
#include <boost/process.hpp>
#include <chrono>
#include <fstream>
#include <sys/resource.h>
#include <thread>
extern "C" {
#include "gen.h"
#include "ralloc.h"
float dryhc08iCode(iCode *ic);
// bool hc08_assignment_optimal;
}

#define REG_A 0
#define REG_X 1
#define REG_H 2

// Code for another ic is generated when generating this one. Mark the other as
// generated.
static void extra_ic_generated(iCode *ic) {
  if (ic->op == '>' || ic->op == '<' || ic->op == LE_OP || ic->op == GE_OP ||
      ic->op == EQ_OP || ic->op == NE_OP || ic->op == '^' || ic->op == '|' ||
      ic->op == BITWISEAND) {
    iCode *ifx;
    if (ifx = ifxForOp(IC_RESULT(ic), ic)) {
      OP_SYMBOL(IC_RESULT(ic))->for_newralloc = false;
      OP_SYMBOL(IC_RESULT(ic))->regType = REG_CND;
      ifx->generated = true;
    }
  }
  if (ic->op == '-' && IS_VALOP(IC_RIGHT(ic)) &&
      operandLitValue(IC_RIGHT(ic)) == 1 &&
      getSize(operandType(IC_RESULT(ic))) == 1 &&
      !isOperandInFarSpace(IC_RESULT(ic)) &&
      isOperandEqual(IC_RESULT(ic), IC_LEFT(ic))) {
    iCode *ifx;
    if (ifx = ifxForOp(IC_RESULT(ic), ic)) {
      OP_SYMBOL(IC_RESULT(ic))->for_newralloc = false;
      OP_SYMBOL(IC_RESULT(ic))->regType = REG_CND;
      ifx->generated = true;
    }
  }
  if (ic->op == GET_VALUE_AT_ADDRESS) {
    iCode *inc;
    if (inc = hasInchc08(IC_LEFT(ic), ic, getSize(operandType(IC_RIGHT(ic)))))
      inc->generated = true;
  }
}

template <class I_t>
static void add_operand_conflicts_in_node(const cfg_node &n, I_t &I) {
  const iCode *ic = n.ic;

  const operand *result = IC_RESULT(ic);
  const operand *left = IC_LEFT(ic);
  const operand *right = IC_RIGHT(ic);

  if (!result || !IS_SYMOP(result))
    return;

  // Todo: Identify more operations that code generation can always handle and
  // exclude them (as done for the z80-like ports).
  if (ic->op == '=')
    return;

  operand_map_t::const_iterator oir, oir_end, oirs;
  boost::tie(oir, oir_end) =
      n.operands.equal_range(OP_SYMBOL_CONST(result)->key);
  if (oir == oir_end)
    return;

  operand_map_t::const_iterator oio, oio_end;

  if (left && IS_SYMOP(left))
    for (boost::tie(oio, oio_end) =
             n.operands.equal_range(OP_SYMBOL_CONST(left)->key);
         oio != oio_end; ++oio)
      for (oirs = oir; oirs != oir_end; ++oirs) {
        var_t rvar = oirs->second;
        var_t ovar = oio->second;
        if (I[rvar].byte < I[ovar].byte)
          boost::add_edge(rvar, ovar, I);
      }

  if (right && IS_SYMOP(right))
    for (boost::tie(oio, oio_end) =
             n.operands.equal_range(OP_SYMBOL_CONST(right)->key);
         oio != oio_end; ++oio)
      for (oirs = oir; oirs != oir_end; ++oirs) {
        var_t rvar = oirs->second;
        var_t ovar = oio->second;
        if (I[rvar].byte < I[ovar].byte)
          boost::add_edge(rvar, ovar, I);
      }
}

// Return true, iff the operand is placed (partially) in r.
template <class G_t>
static bool operand_in_reg(const operand *o, reg_t r, const i_assignment_t &ia,
                           unsigned short int i, const G_t &G) {
  if (!o || !IS_SYMOP(o))
    return (false);

  if (r >= port->num_regs)
    return (false);

  operand_map_t::const_iterator oi, oi_end;
  for (boost::tie(oi, oi_end) =
           G[i].operands.equal_range(OP_SYMBOL_CONST(o)->key);
       oi != oi_end; ++oi)
    if (oi->second == ia.registers[r][1] || oi->second == ia.registers[r][0])
      return (true);

  return (false);
}

template <class G_t, class I_t>
static bool operand_is_ax(const operand *o, const assignment &a,
                          unsigned short int i, const G_t &G, const I_t &I) {
  if (!o || !IS_SYMOP(o))
    return (false);

  operand_map_t::const_iterator oi, oi2, oi_end;
  boost::tie(oi, oi_end) = G[i].operands.equal_range(OP_SYMBOL_CONST(o)->key);

  if (oi == oi_end)
    return (false);

  oi2 = oi;
  oi2++;
  if (oi2 == oi_end)
    return (false);

  // Register combinations code generation cannot handle yet (AX, AH, XH, HA).
  if (std::binary_search(a.local.begin(), a.local.end(), oi->second) &&
      std::binary_search(a.local.begin(), a.local.end(), oi2->second)) {
    const reg_t l = a.global[oi->second];
    const reg_t h = a.global[oi2->second];
    if (l == REG_X && h == REG_A)
      return (true);
  }

  return (false);
}

static void write_into_csv(float c, int i, int time) {
  std::ofstream outputFile("optimalCost.csv", std::ios_base::app);
  if (outputFile.is_open()) { // Check if the file was successfully opened
    // Write some text into the file
    outputFile << std::string(dstFileName) << "," << currFunc->name << "," << c
               << "," << time << "," << i
               << "\n"; // Write a line of text to the file
    // Close the file
    outputFile.close(); // Close the file after writing

    std::cout << "Text has been written to the file."
              << std::endl; // Display a success message
  } else {
    std::cout << "Failed to create the file."
              << std::endl; // Display an error message if file creation failed
  }
}

static void write_into_csv(int i, int time) {
  std::ofstream outputFile("optimalc_time.csv", std::ios_base::app);
  if (outputFile.is_open()) { // Check if the file was successfully opened
    // Write some text into the file
    outputFile << std::string(dstFileName) << "," << time << "," << i
               << "\n"; // Write a line of text to the file
    // Close the file
    outputFile.close(); // Close the file after writing

    std::cout << "Text has been written to the file."
              << std::endl; // Display a success message
  } else {
    std::cout << "Failed to create the file."
              << std::endl; // Display an error message if file creation failed
  }
}

struct Coloring {
  static constexpr char MAX_REGS = 30;
  // optimize this?
  boost::container::small_vector<std::pair<var_t, char>, 4> coloring;

  bool operator==(const Coloring &other) const {
    return coloring == other.coloring;
  }

  bool operator<(const Coloring &other) const {
    return coloring < other.coloring;
  }

  void canonicalize() {
    char reg_map[MAX_REGS];
    std::fill_n(reg_map, MAX_REGS, -1);

    std::sort(coloring.begin(), coloring.end(),
              [](const std::pair<var_t, int> &a,
                 const std::pair<var_t, int> &b) { return a.first < b.first; });

    char current_reg = 0;
    for (auto &p : coloring) {
      if (reg_map[p.second] == -1)
        reg_map[p.second] = current_reg++;
      p.second = reg_map[p.second];
    }
  }

  // true if successful
  bool merge(const Coloring &other, Coloring &result) const {
    char reg_map[MAX_REGS];
    std::fill_n(reg_map, MAX_REGS, -1);

    bool used[MAX_REGS] = {0};
    decltype(coloring) tmp;
    for (const auto &p : other.coloring) {
      auto matching =
          std::lower_bound(coloring.cbegin(), coloring.cend(), p.first,
                           [](const std::pair<var_t, char> &elem, var_t value) {
                             return elem.first < value;
                           });
      // we don't really care if it does not appear in both maps
      if (matching == coloring.cend()) {
        tmp.push_back(p);
        continue;
      }

      if (reg_map[p.second] == -1) {
        // we require the mapping to be a permutation
        // std::cout << static_cast<int>(p.second) << "->" <<
        // static_cast<int>(matching->second) << std::endl;
        if (used[matching->second])
          return false;
        reg_map[p.second] = matching->second;
        used[matching->second] = true;
      } else if (reg_map[p.second] != matching->second) {
        return false;
      }
    }

    for (auto &p : tmp) {
      char r = reg_map[p.second];
      if (r == -1) {
        while (used[++r] && r < port->num_regs)
          continue;

        // running out of registers...
        if (used[r])
          return false;
        used[r] = true;
        reg_map[p.second] = r;
      }
      p.second = r;
    }

    result.coloring.resize(coloring.size() + tmp.size());
    std::merge(
        coloring.begin(), coloring.end(), tmp.begin(), tmp.end(),
        result.coloring.begin(),
        [](const std::pair<var_t, char> &a, const std::pair<var_t, char> &b) {
          return a.first < b.first;
        });
    result.canonicalize();
    return true;
  }
};

template <> struct std::hash<Coloring> {
  std::size_t operator()(const Coloring &c) const {
    return boost::hash_value(c.coloring);
  }
};

namespace {
// reuse the dom tree part
class DomTree {
public:
  DomTree() = default;

  template <class G_t> DomTree(const G_t &G, unsigned int start = 0) {
    GetPostorder(G, start);
    GetVertToRevPostorder();
    GetIDom(G, start);
  }

  template <class G_t>
  std::set<unsigned int> GetUnreachable(const G_t &G) const {
    std::set<unsigned int> unreachable;
    for (unsigned int i = 0; i < boost::num_vertices(G); i++) {
      if (std::find(postorder_.begin(), postorder_.end(), i) ==
          postorder_.end()) {
        unreachable.insert(i);
      }
    }
    return unreachable;
  }

  bool Dominates(unsigned int u, unsigned int v) const {
    v = idom_[v];
    unsigned int rpo1 = vertToRpo_[u];
    unsigned int rpo2 = vertToRpo_[v];

    while (rpo1 < rpo2) {
      v = idom_[v];
      rpo2 = vertToRpo_[v];
    }
    return rpo1 == rpo2;
  }

  unsigned int GetIDom(unsigned int u) const { return idom_[u]; }

private:
  std::vector<unsigned int> idom_;
  std::vector<unsigned int> postorder_;
  std::vector<unsigned int> vertToRpo_;

  template <class G_t> void GetPostorder(const G_t &G, unsigned int start) {
    typedef
        typename boost::graph_traits<G_t>::adjacency_iterator adjacency_iter_t;
    struct State {
      unsigned int id;
      adjacency_iter_t curr;
      adjacency_iter_t end;
    };
    auto index = boost::get(boost::vertex_index, G);

    std::vector<bool> visited(boost::num_vertices(G), false);
    std::vector<State> stack;

    auto adj_0 = boost::adjacent_vertices(start, G);
    visited[start] = true;
    stack.push_back({start, adj_0.first, adj_0.second});

    while (!stack.empty()) {
      State v = stack.back();
      bool found_unvisited = false;
      for (; v.curr != v.end; ++v.curr) {
        unsigned int u = index[*v.curr];
        if (!visited[u]) {
          visited[u] = true;
          auto adj_u = boost::adjacent_vertices(u, G);
          stack.push_back({u, adj_u.first, adj_u.second});
          found_unvisited = true;
          break;
        }
      }
      if (!found_unvisited) {
        postorder_.push_back(v.id);
        stack.pop_back();
      }
    }
  }

  void GetVertToRevPostorder() {
    vertToRpo_.resize(postorder_.size(), INT_MAX);
    for (unsigned int i = 0; i < postorder_.size(); i++)
      vertToRpo_[postorder_[i]] = postorder_.size() - i - 1;
  }

  // This is an implementation of the algorithm described in
  //
  //   A Simple, Fast Dominance Algorithm
  //   Keith D. Cooper, Timothy J. Harvey, and Ken Kennedy
  //   Department of Computer Science, Rice University, Houston, Texas, USA
  //   TR-06-33870
  //   https://www.cs.rice.edu/~keith/EMBED/dom.pdf
  template <class G_t> void GetIDom(const G_t &G, unsigned start) {
    auto index = boost::get(boost::vertex_index, G);
    auto &vertToRpo = vertToRpo_;
    auto &idom = idom_;
    idom.resize(boost::num_vertices(G), INT_MAX);
    auto intersect = [&vertToRpo, &idom](unsigned int a, unsigned int b) {
      unsigned int finger1 = a;
      unsigned int finger2 = b;
      while (finger1 != finger2) {
        unsigned int rpo1 = vertToRpo[finger1];
        unsigned int rpo2 = vertToRpo[finger2];
        if (rpo1 > rpo2)
          finger1 = idom[finger1];
        else
          finger2 = idom[finger2];
      }
      return finger1;
    };

    idom[start] = start;

    auto inv_graph = boost::make_reverse_graph(G);
    bool changed = true;
    while (changed) {
      changed = false;
      for (int i = postorder_.size() - 1; i >= 0; i--) {
        unsigned int v = postorder_[i];
        // skip starting vertex
        if (v == start)
          continue;
        unsigned int new_idom = INT_MAX;
        // predecessors
        for (auto it : boost::make_iterator_range(
                 boost::adjacent_vertices(v, inv_graph))) {
          unsigned int u = index[it];
          if (idom[u] == INT_MAX)
            continue;
          if (new_idom == INT_MAX)
            new_idom = u;
          else
            new_idom = intersect(new_idom, u);
        }
        if (idom[v] != new_idom) {
          idom[v] = new_idom;
          changed = true;
        }
      }
    }
  }
};

std::vector<cfg_node> cfg_nodes;

struct SPL_tree_node {
  // index to each special node in our decomposition
  int begin;
  int end;

  boost::container::small_flat_set<int, 4>
      alive_variables; // all alive variable at special node

  int type; // to show the type of this node, 0 is atomic, 1 is series, 2 is
            // parallel, 3 is loop

  int parent = -1;
  bool ifleft = true;
  int left;
  int right;
  // bool in_loop = false;

  SPL_tree_node(int b, int e) {
    begin = b;
    end = e;
    for (auto v : cfg_nodes[b].alive) {
      alive_variables.insert(v);
    }

    for (auto v : cfg_nodes[e].alive) {
      alive_variables.insert(v);
    }
  }
};

std::vector<SPL_tree_node>
    SPL_tree_nodes; // a vector contains all SPL tree nodes

// add the SPLnode to the SPL_tree_nodes
void update_list(SPL_tree_node current) {
  // check if it is the root
  if (SPL_tree_nodes.empty()) {
    SPL_tree_nodes.push_back(current);
    return;
  }

  // update the parent node
  if (current.ifleft) {
    SPL_tree_nodes[current.parent].left = SPL_tree_nodes.size();
  } else {
    SPL_tree_nodes[current.parent].right = SPL_tree_nodes.size();
  }

  SPL_tree_nodes.push_back(current);
}

// 1. I reused the unreachable node part and remove all the unreachable nodes
// 2. Get the backedge with domtree, and generate all loop node directly, store
// in loops
// 3. Remove the back edges and the exit edges from the loop
// 4. Go through the CFG and find all parallel part, generate the parallel nodes
// 5. go through the CFG again and convert the cfg to spl_decomposition
void convert_cfg_to_spl(cfg_t &cfg) {
  SPL_tree_nodes.clear();
  cfg_nodes.clear();

  for (int v = 0; v < boost::num_vertices(cfg); v++) {
    cfg_nodes.push_back(cfg[v]);
  }
  cfg_t simplified_;
  std::set<unsigned int> unreachable_;

  boost::copy_graph(cfg, simplified_);
  auto index = boost::get(boost::vertex_index, simplified_);

  DomTree domtree(simplified_);

  // remove all unreachable nodes
  int lastVertex = boost::num_vertices(simplified_) - 1;
  unreachable_ = domtree.GetUnreachable(simplified_);
  if (!unreachable_.empty()) {
    for (unsigned int v : unreachable_) {
      boost::clear_out_edges(v, simplified_);
    }
    // recalculate dominator tree
    domtree = DomTree(simplified_);
  }
  while (unreachable_.find(lastVertex) != unreachable_.end()) {
    lastVertex--;
  }
  // get and store all loop and related special nodes

  auto backedges = boost::container::flat_map<
      unsigned int, boost::container::small_vector<unsigned int, 20>>();

  for (unsigned int v = 0; v <= lastVertex; v++) {
    for (auto iter : boost::make_iterator_range(
             boost::inv_adjacent_vertices(v, simplified_))) {
      unsigned int u = index[iter];
      // edge: u -> v
      // if v dominates u, u -> v is a backedge
      if (domtree.Dominates(v, u)) {
        auto entry = backedges.find(v);
        if (entry == backedges.end())
          backedges.insert(
              std::pair<unsigned int,
                        boost::container::small_vector<unsigned int, 20>>(v,
                                                                          {u}));
        else
          (*entry).second.push_back(u);
      }
    }
  }
  auto loops = boost::container::flat_map<unsigned int, SPL_tree_node>();
  for (auto iter : backedges) {
    auto end = std::max_element(iter.second.begin(), iter.second.end());
    SPL_tree_node node(iter.first, *end);
    node.type = 3;
    for (auto con_iter : iter.second) {
      // if (con_iter < *end)
      // {
      //   node.con.push_back(con_iter);
      remove_edge(con_iter, iter.first, simplified_);
      // std::cout << "continue edge remove: " << con_iter << "->" << iter.first
      // << std::endl;
      //}
    }

    auto breaks = std::vector<int>();

    for (int v = iter.first; v <= *end; v++) {
      for (auto iterb : boost::make_iterator_range(
               boost::adjacent_vertices(v, simplified_))) {
        if (index[iterb] > *end) {
          //   node.brk.push_back(index[iterb]);
          remove_edge(v, index[iterb], simplified_);
          // check if v in breaks
          breaks.push_back(index[iterb]);

          // std::cout << "break edge remove: " << v << "->" << index[iterb] <<
          // std::endl;
          break;
        }
      }
    }
    // find max v in breaks
    if (!breaks.empty()) {
      auto min_v = std::min_element(breaks.begin(), breaks.end());
      boost::add_edge(*end, *min_v, simplified_);
    }

    loops.insert(std::pair<unsigned int, SPL_tree_node>(iter.first, node));
  }

  // get all parallel part
  auto parallels = boost::container::flat_map<unsigned int, SPL_tree_node>();
  for (unsigned int v = 0; v <= lastVertex; v++) {
    if (boost::out_degree(v, simplified_) > 1) {
      int p_num = 1;
      int current = v;
      while (p_num != 0 && v <= lastVertex) {
        auto adj_v = boost::adjacent_vertices(current, simplified_);
        if (adj_v.first == adj_v.second) {
          break;
        }
        int next = index[*(adj_v.first)];
        if (boost::in_degree(next, simplified_) == 2) {
          p_num--;
        }
        if (boost::out_degree(next, simplified_) == 2) {
          p_num++;
        }
        current = next;
      }
      if (p_num != 0) {
        throw std::runtime_error("invalid cfg due to p_num check");
      }
      SPL_tree_node node(v, current);
      node.type = 2;
      parallels.insert(std::pair<unsigned int, SPL_tree_node>(v, node));
    }
  }

  // all nodes that need to be analysised and separated are in the worklist
  std::vector<SPL_tree_node> worklist;
  SPL_tree_node root(0, lastVertex);
  worklist.push_back(root);
  while (!worklist.empty()) {
    SPL_tree_node current = worklist.back();
    worklist.pop_back();
    if (current.begin > current.end) {
      throw std::runtime_error(
          "invalid cfg due to current.begin > current.end");
    }
    // if begin and end is the same, then it is atomic
    if (current.begin == current.end) {
      current.type = 0;
      update_list(current);
      continue;
    }

    auto adj_v_range = boost::adjacent_vertices(current.begin, simplified_);
    if (adj_v_range.first == adj_v_range.second) {
      throw std::runtime_error("invalid cfg due to no adjacent vertices");
    }

    // if the begin is the begin of a loop.
    if (loops.find(current.begin) != loops.end()) {

      SPL_tree_node loop_node = loops.find(current.begin)->second;
      loops.erase(current.begin);
      if (loop_node.end != current.end) {
        auto next_loop = boost::adjacent_vertices(loop_node.end, simplified_);
        if (next_loop.first == next_loop.second) {
          throw std::runtime_error("invalid cfg due to no next loop");
        }
        SPL_tree_node after_loop(index[*(next_loop.first)], current.end);
        after_loop.parent = SPL_tree_nodes.size();
        after_loop.ifleft = false;
        worklist.push_back(after_loop);
        loop_node.parent = SPL_tree_nodes.size();
        current.type = 1;
        update_list(current);
      } else {

        loop_node.parent = current.parent;
      }

      update_list(loop_node);
      SPL_tree_node loop_body(loop_node.begin, loop_node.end);
      loop_body.parent = SPL_tree_nodes.size() - 1;
      worklist.push_back(loop_body);

      continue;
    }

    // if the begin vertex begin a parallel
    if (parallels.find(current.begin) != parallels.end()) {

      SPL_tree_node parallel_node = parallels.find(current.begin)->second;
      parallels.erase(current.begin);
      if (parallel_node.end != current.end) {
        auto next_parallel =
            boost::adjacent_vertices(parallel_node.end, simplified_);
        if (next_parallel.first == next_parallel.second) {
          throw std::runtime_error("invalid cfg due to no next parallel");
        }
        SPL_tree_node after_parallel(index[*(next_parallel.first)],
                                     current.end);
        after_parallel.parent = SPL_tree_nodes.size();
        after_parallel.ifleft = false;
        worklist.push_back(after_parallel);
        parallel_node.parent = SPL_tree_nodes.size();
        current.type = 1;
        update_list(current);
      } else {
        parallel_node.parent = current.parent;
      }
      update_list(parallel_node);

      std::vector<int> adj_v;
      for (auto iter : boost::make_iterator_range(adj_v_range)) {
        adj_v.push_back(index[iter]);
      }
      if (adj_v.size() != 2) {
        throw std::runtime_error(
            "invalid CFG due to wrong degree of parallel node");
      }

      SPL_tree_node left(current.begin, parallel_node.end);
      SPL_tree_node right(parallel_node.begin, parallel_node.end);
      int head = current.begin;
      left.parent = SPL_tree_nodes.size() - 1;
      right.parent = SPL_tree_nodes.size() - 1;
      right.ifleft = false;
      left.type = 1;
      right.type = 1;

      update_list(left);
      SPL_tree_node left_head(head, head);
      left_head.parent = SPL_tree_nodes.size() - 1;
      left_head.type = 0;
      SPL_tree_node left_body(adj_v[0], parallel_node.end);
      left_body.parent = SPL_tree_nodes.size() - 1;
      left_body.ifleft = false;
      update_list(left_head);
      worklist.push_back(left_body);

      update_list(right);
      SPL_tree_node right_head(head, head);
      right_head.parent = SPL_tree_nodes.size() - 1;
      right_head.type = 0;
      SPL_tree_node right_body(adj_v[1], parallel_node.end);
      right_body.parent = SPL_tree_nodes.size() - 1;
      update_list(right_head);
      right_body.ifleft = false;
      worklist.push_back(right_body);
      continue;
    }

    // it is simple series type

    current.type = 1;

    SPL_tree_node head_part(current.begin, current.begin);

    SPL_tree_node right_part(index[*(adj_v_range.first)], current.end);
    update_list(current);

    head_part.parent = SPL_tree_nodes.size() - 1;
    head_part.type = 0;
    right_part.parent = SPL_tree_nodes.size() - 1;
    right_part.ifleft = false;
    update_list(head_part);

    worklist.push_back(right_part);
  }
}

using Entry = std::set<Coloring>;

template <typename T>
void handle_atom(Entry &assignments, const T &live_vars,
                 const con_t &conflict_graph) {
  Coloring c;
  for (auto x : live_vars)
    c.coloring.push_back(std::make_pair(x, 0));
  c.canonicalize();

  if (live_vars.empty()) {
    assignments.insert(c);
    return;
  }

  std::function<void(int, int)> rec = [&](int i, int max_r) {
    for (int r = 0; r <= std::min(max_r, port->num_regs - 1); r++) {
      bool conflict = false;
      for (int j = 0; j < i; j++)
        if (c.coloring[j].second == r &&
            boost::edge(c.coloring[i].first, c.coloring[j].first,
                        conflict_graph)
                .second) {
          conflict = true;
          break;
        }
      if (conflict)
        continue;
      c.coloring[i].second = r;
      if (i + 1 == c.coloring.size())
        assignments.insert(c);
      else
        rec(i + 1, std::max(r + 1, max_r));
    }
  };
  rec(0, 0);
}

template <typename T>
void handle_composition(Entry &result, const Entry &lhs, const Entry &rhs,
                        const T &live_vars) {
  Coloring merged, final;
  for (const auto &f1 : lhs) {
    for (const auto &f2 : rhs) {
      if (f1.merge(f2, merged)) {
        final.coloring.clear();
        for (const auto &p : merged.coloring) {
          if (live_vars.find(p.first) != live_vars.end())
            final.coloring.push_back(p);
        }
        result.insert(final);
      }
    }
  }
}

// just restriction...
template <typename T>
void handle_loop(Entry &result, const Entry &body, const T &live_vars) {
  Coloring final;
  for (const auto &f : body) {
    final.coloring.clear();
    for (const auto &p : f.coloring) {
      if (live_vars.find(p.first) != live_vars.end())
        final.coloring.push_back(p);
    }
    result.insert(final);
  }
}

int cfg_decomp_regalloc(const cfg_t &control_flow_graph,
                        const con_t &conflict_graph) {
  int old_regs = port->num_regs;
  bool ok = true;
  int r;
  for (r = 0; r < 21; r++) {
    port->num_regs = r;

    ok = true;

    std::vector<Entry> entries(SPL_tree_nodes.size());
    for (int i = SPL_tree_nodes.size() - 1; i >= 0; i--) {
      auto &entry = entries[i];
      const auto &node = SPL_tree_nodes[i];
      const auto &live = node.alive_variables;
      switch (node.type) {
      case 0:
        handle_atom(entry, live, conflict_graph);
        break;
      case 1:
      case 2:
        handle_composition(entry, entries[node.left], entries[node.right],
                           live);
        entries[node.left].clear();
        entries[node.right].clear();
        break;
      case 3:
        handle_loop(entry, entries[node.left], live);
        entries[node.left].clear();
        break;
      }
      if (entry.empty()) {
        ok = false;
        break;
      }
    }
    if (ok) {
      break;
    }
  }
  port->num_regs = old_regs;
  if (!ok) {
    return -1;
  }

  return r;
}
} // namespace

float hc08_ralloc3_cc(ebbIndex *ebbi) {

#ifdef DEBUG_RALLOC_DEC
  std::cout << "Processing " << currFunc->name << " from " << dstFileName
            << "\n";
  std::cout.flush();
#endif
  char * IN_BUILD = getenv("IN_BUILD");
  if (IN_BUILD != NULL && std::string(IN_BUILD) == "1")
    return -1;

  cfg_t control_flow_graph;

  con_t conflict_graph;

  iCode *ic = create_cfg(control_flow_graph, conflict_graph, ebbi);

  // ------------------------------------------------------------------------------------------------
  struct rlimit limits;
  // 4GB limit
  limits.rlim_cur = 16l * 1024 * 1024 * 1024;
  limits.rlim_max = limits.rlim_cur;
  setrlimit(RLIMIT_AS, &limits);
  // ------------------------------------------------------------------------------------------------

  constexpr int max_w = 21;

  auto benchmark = [&](std::function<void()> f) {
    auto start = std::chrono::high_resolution_clock::now();
    f();
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(end - start)
        .count();
  };

  // linear scan
  int linear_scan_r = -1;
  auto linear_scan_time = benchmark([&]() {
    for (int i = 0; i < boost::num_vertices(control_flow_graph); i++) {
      linear_scan_r =
          std::max(linear_scan_r, (int)control_flow_graph[i].alive.size());
    }
  });

  // graph coloring with minisat
  int graph_color_r = -1;
  auto graph_color_time = benchmark([&]() {
    if (linear_scan_r == 0) {
      graph_color_r = 0;
      return;
    }
    for (int w = 1; w <= max_w; w++) {
      namespace bp = boost::process;
      bp::ipstream out;
      int num_clauses = 0;
      int n = boost::num_vertices(conflict_graph);

      std::ofstream minisat_in;
      minisat_in.open("./tmp.txt");

      int edges = 0;
      for (int i = 0; i < n; i++) {
        for (auto it : boost::make_iterator_range(
                 boost::adjacent_vertices(i, conflict_graph))) {
          edges++;
        }
      }

      minisat_in << "p cnf " << (n * w) << " "
                 << (n * w * (w - 1) / 2 + n + w * edges) << std::endl;
      // for each variable, we have w variables representing w registers
      // satisfying two requirements
      // 1. at least one should be true
      // 2. no pair should be true simultaneously
      // so there are w*(w-1)/2 + 1 clauses per variable to encode single color
      // per variable
      for (int i = 0; i < n; i++) {
        // at least one
        for (int j = 0; j < w; j++)
          minisat_in << (i * w + j + 1) << " ";
        minisat_in << "0" << std::endl;

        // at most one
        for (int j = 1; j < w; j++)
          for (int k = 0; k < j; k++)
            minisat_in << -(i * w + j + 1) << " " << -(i * w + k + 1) << " 0"
                       << std::endl;
      }

      // for each pair of conflicting variables, we have w clauses,
      // one for each color
      for (int i = 0; i < n; i++) {
        for (auto j : boost::make_iterator_range(
                 boost::adjacent_vertices(i, conflict_graph))) {
          for (int k = 0; k < w; k++)
            minisat_in << -(i * w + k + 1) << " " << -(int(j) * w + k + 1)
                       << " 0" << std::endl;
        }
      }

      minisat_in.close();

      // tried using stdin, somehow opstream doesn't work
      bp::child c("minisat -verb=0 -mem-lim=4096 -cpu-lim=10 ./tmp.txt",
                  bp::std_out > out);

      c.wait();
      std::string value(std::istreambuf_iterator<char>(out), {});
      if (value.find("UNSATISFIABLE") == -1) {
        // not satisfiable either...
        if (value.find("SATISFIABLE") == -1) {
          graph_color_r = -1;
          break;
        }
        graph_color_r = w;
        break;
      }
    }
  });

  int spl_r = -1;
  auto spl_time = benchmark([&]() {
    try {
      convert_cfg_to_spl(control_flow_graph);
      spl_r = cfg_decomp_regalloc(control_flow_graph, conflict_graph);
    } catch (std::runtime_error &err) {
      std::cout << err.what() << std::endl;
      spl_r = -1;
    }
  });

  int treedec_r = -1;
  auto treedec_time = benchmark([&]() {
    tree_dec_t tree_decomposition;

    get_nice_tree_decomposition(tree_decomposition, control_flow_graph);

    alive_tree_dec(tree_decomposition, control_flow_graph);

    good_re_root(tree_decomposition);
    nicify(tree_decomposition);
    alive_tree_dec(tree_decomposition, control_flow_graph);
    tree_dec_var_t intersection_treedec =
        intersection_treedecomp(control_flow_graph, tree_decomposition);
    auto root = find_root(tree_decomposition);

    int w = 0;
    for (; w <= 10; w++) {
      try {
        phmap::parallel_flat_hash_set<coloring_t> colorings;
        restart_timer();
        if (try_coloring_treedec(w, root, intersection_treedec, conflict_graph,
                                 colorings)) {
          treedec_r = w;
          break;
        }
      } catch (const std::runtime_error &e) {
        if (strcmp(e.what(), "timeout") == 0) {
          std::cout << "timeout doing tree decomposition with w = " << w
                    << std::endl;
        }
        return;
      } catch (const std::bad_alloc &e) {
        std::cout << "OOM when w = " << w << std::endl;
        return;
      }
    }
  });

  int pathdec_r = -1;
  auto pathdec_time = benchmark([&]() {
    auto decomp = intersection_decomp(control_flow_graph);
    int w = 0;
    for (; w <= 10; w++) {
      try {
        restart_timer();
        if (try_coloring(w, decomp, conflict_graph)) {
          pathdec_r = w;
          break;
        }
      } catch (const std::runtime_error &e) {
        if (strcmp(e.what(), "timeout") == 0) {
          std::cout << "timeout doing path decomposition with w = " << w
                    << std::endl;
        }
        return;
      } catch (const std::bad_alloc &e) {
        std::cout << "OOM when w = " << 2 << std::endl;
        return;
      }
    }
  });

  int puzzle_r = -1;
  auto puzzle_time = benchmark([&]() {
    int w = 0;
    for (; w <= max_w; w++) {
      try {
        restart_timer();
        if (getGlobalReg(w, control_flow_graph)) {
          puzzle_r = w;
          break;
        }
      } catch (const std::runtime_error &e) {
        if (strcmp(e.what(), "timeout") == 0) {
          std::cout << "timeout solving puzzle with w = " << w << std::endl;
        }
        return;
      }
    }
  });

  std::cout << "linear_scan_r: " << linear_scan_r << ", "
            << "graph_color_r: " << graph_color_r << ", "
            << "spl_r: " << spl_r << ", "
            << "treedec_r: " << treedec_r << ", "
            << "pathdec_r: " << pathdec_r << ", "
            << "puzzle_r: " << puzzle_r << ", "
            << "linear_scan_time: " << linear_scan_time << ", "
            << "graph_color_time: " << graph_color_time << ", "
            << "spl_time: " << spl_time << ", "
            << "treedec_time: " << treedec_time << ", "
            << "pathdec_time: " << pathdec_time << ", "
            << "puzzle_time: " << puzzle_time
            << std::endl;

  return -1;
}
