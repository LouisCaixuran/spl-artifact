#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/transpose_graph.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>

// much faster than std::set due to less allocation cost
typedef boost::container::flat_set<
    unsigned int, std::less<int>,
    boost::container::small_vector<unsigned int, 20>>
    PathDecBag;

// Verify that the path decomposition is a valid one
// Note that this does not check if it is a nice path decomposition, but that
// should be trivial to check
template <class G_t>
bool VerifyPathDecomposition(const G_t &G, std::vector<PathDecBag> &dec,
                             std::set<unsigned int> &unreachable) {
  auto index = boost::get(boost::vertex_index, G);
  std::vector<std::pair<unsigned int, unsigned int>> vertexRange(
      boost::num_vertices(G), {UINT32_MAX, UINT32_MAX});
  for (unsigned int j = 0; j < dec.size(); j++) {
    for (unsigned int i : dec[j]) {
      if (vertexRange[i].first == UINT32_MAX) {
        vertexRange[i].first = j;
        vertexRange[i].second = j;
      } else if (vertexRange[i].second == j - 1) {
        vertexRange[i].second = j;
      } else {
        std::cout << "disconnected interval" << std::endl;
        return false;
      }
    }
  }

  for (unsigned int i = 0; i < boost::num_vertices(G); i++) {
    // skip unreachable nodes
    if (unreachable.find(i) != unreachable.end())
      continue;
    // assert the vertex must belong to at least one bag
    if (vertexRange[i].first == UINT32_MAX) {
      std::cout << "vertex " << i << " does not belong to any bag" << std::endl;
      return false;
    }
    for (auto j_curr :
         boost::make_iterator_range(boost::adjacent_vertices(i, G))) {
      unsigned int j = index[j_curr];
      // j must occur in one of the bags in the interval
      bool j_occured = false;
      for (unsigned int k = vertexRange[i].first; k <= vertexRange[i].second;
           k++) {
        if (std::find(dec[k].begin(), dec[k].end(), j) != dec[k].end()) {
          j_occured = true;
          break;
        }
      }
      if (!j_occured) {
        std::cout << "edge " << i << "->" << j << " does not belong to any bag"
                  << std::endl;
        return false;
      }
    }
  }
  return true;
}

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

template <class G_t> class Cfg {
public:
  Cfg(const G_t &G) : handled_(boost::num_vertices(G), false) {
    boost::copy_graph(G, simplified_);
    RemoveBreakContinue();

    auto index = boost::get(boost::vertex_index, simplified_);
    bagStore_.reserve(2 * (lastVertex + 1));
    stack_.reserve(20);
    AddTask(0, lastVertex);
    PathDecStep(index);
    ApplyPatches(0);

    // Make the vertex range continuous.
    //
    // Problem: we cannot handle nodes with out-degree > 2 very well
    // Consider v -> u -> a -> b, v -> w -> a -> b, v -> x -> b
    // The immediate postdominator of v is b, but the immediate postdominator of
    // u and w is a. Our algorithm will not know that a needs to be duplicated
    // and will produce an incorrect path decomposition. We can handle it
    // appropriately in the algorithm but that would be slow, and this way of
    // handling it will not cause a worse path decomposition so it should be
    // fine.
    std::vector<std::pair<unsigned int, unsigned int>> vertexRange(
        lastVertex + 1, {UINT32_MAX, UINT32_MAX});
    for (unsigned int j = 0; j < path_.size(); j++) {
      for (unsigned int i : bagStore_[path_[j]]) {
        if (vertexRange[i].first == UINT32_MAX) {
          vertexRange[i].first = j;
          vertexRange[i].second = j;
        } else if (vertexRange[i].second == j - 1) {
          vertexRange[i].second = j;
        } else {
          // fix discontinuous
          for (unsigned int k = vertexRange[i].second; k < j; k++)
            bagStore_[path_[k]].insert(i);
          vertexRange[i].second = j;
        }
      }
    }

    // Add back the removed edges
    for (auto iter : fixup_) {
      unsigned int u = iter.first;
      // find last occurance
      auto first = vertexRange[u].first;
      auto last = vertexRange[u].second;
      for (auto v : iter.second) {
        auto vFirst = vertexRange[v].first;
        auto vLast = vertexRange[v].second;
        assert(first != UINT32_MAX);
        assert(vFirst != UINT32_MAX);
        if (last < vFirst) {
          // u is before v
          for (unsigned int i = last; i < vFirst; i++)
            bagStore_[path_[i]].insert(v);
          vertexRange[v].first = last;
        } else if (first > vLast) {
          // u is after v
          for (unsigned int i = vLast + 1; i <= first; i++)
            bagStore_[path_[i]].insert(v);
          vertexRange[v].second = first;
        }
      }
    }
    nicify();
  }

  std::vector<PathDecBag> &GetPathDec() { return nicePath_; }

  std::set<unsigned int> &GetUnreachable() { return unreachable_; }

private:
  typedef
      typename boost::graph_traits<G_t>::adjacency_iterator adjacency_iter_t;

  // a struct that stores intermediate states of path decomposition, used to
  // manually maintain a stack of states instead of doing recursion on the graph
  struct PathDecState {
    adjacency_iter_t iter_start;
    adjacency_iter_t iter_curr;
    adjacency_iter_t iter_end;
    unsigned int start;
    unsigned int end;
    unsigned int iend;
    unsigned int index1;
    unsigned int index2;
    bool swapped;
  };

  G_t simplified_;
  DomTree postdomtree_;
  boost::container::flat_map<unsigned int, std::set<unsigned int>> fixup_;
  unsigned int addOnce = UINT32_MAX;
  std::vector<bool> handled_;
  std::vector<PathDecBag> bagStore_;
  std::vector<unsigned int> path_;
  std::vector<PathDecBag> nicePath_;
  std::vector<PathDecState> stack_;
  std::set<unsigned int> unreachable_;

  std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> patches_;
  unsigned int lastVertex = UINT32_MAX;

  // TODO: optimize this part, perhaps we have to change the API altogether
  void nicify() {
    nicePath_ = std::vector<PathDecBag>();
    nicePath_.reserve(2 * (lastVertex + 1));
    nicePath_.push_back({});
    for (auto bag_index : path_) {
      auto &bag = bagStore_[bag_index];
      boost::container::small_vector<unsigned int, 20> tmp;
      std::set_difference(nicePath_.back().begin(), nicePath_.back().end(),
                          bag.begin(), bag.end(), std::back_inserter(tmp));
      auto newBag = PathDecBag(nicePath_.back());
      for (auto v : tmp) {
        newBag.erase(v);
        nicePath_.emplace_back(newBag);
      }
      tmp.clear();
      std::set_difference(bag.begin(), bag.end(), nicePath_.back().begin(),
                          nicePath_.back().end(), std::back_inserter(tmp));
      for (auto v : tmp) {
        newBag.insert(v);
        nicePath_.emplace_back(newBag);
      }
    }
    while (!nicePath_.back().empty()) {
      auto newBag = PathDecBag(nicePath_.back());
      newBag.erase(*newBag.begin());
      nicePath_.emplace_back(newBag);
    }
  }

  template <class I_t> void PathDecStep(const I_t &index) {
    unsigned int lastSubsumption = UINT32_MAX;
    while (1) {
      PathDecState &state = stack_.back();
      unsigned int &start = state.start;
      if (start == state.end) {
        AddBag<1>({state.start});
        stack_.pop_back();
        if (stack_.empty())
          return;
        continue;
      }

      unsigned int newStart = UINT32_MAX;
      // optimize short-circuiting boolean by subsumption...
      if (state.iter_start == state.iter_curr &&
          std::distance(state.iter_start, state.iter_end) == 2) {
        unsigned int u = index[*state.iter_start];
        unsigned int v = index[*(state.iter_start + 1)];
        // u/v subsumes start
        if (boost::edge(u, v, simplified_).second) {
          AddBag<3>({start, u, v});
          newStart = u;
          lastSubsumption = v;
        } else if (boost::edge(v, u, simplified_).second) {
          AddBag<3>({start, u, v});
          newStart = v;
          lastSubsumption = v;
        }
      }

      if (newStart == UINT32_MAX)
        lastSubsumption = UINT32_MAX;

      // reorder 2 branches
      if (state.iter_curr == state.iter_end) {
        if (std::distance(state.iter_start, state.iter_end) == 2) {
          assert(state.index1 != UINT32_MAX);
          assert(state.index2 != UINT32_MAX);
          ApplyPatches(state.index1);
          unsigned int width1 = 0;
          unsigned int width2 = 0;
          for (unsigned int i = state.index1; i < state.index2; i++)
            width1 = std::max(width1, (unsigned int)bagStore_[path_[i]].size());
          for (unsigned int i = state.index2; i < path_.size(); i++)
            width2 = std::max(width2, (unsigned int)bagStore_[path_[i]].size());

          unsigned int size1 = 0, size2 = 0;
          if (state.index1 > 1 && state.index2 != path_.size()) {
            // check continuous
            boost::container::small_vector<unsigned int, 20> buffer;
            std::set_intersection(bagStore_[path_[state.index1 - 1]].begin(),
                                  bagStore_[path_[state.index1 - 1]].end(),
                                  bagStore_[path_[state.index1]].begin(),
                                  bagStore_[path_[state.index1]].end(),
                                  std::back_inserter(buffer));
            size1 = buffer.size();
            if (std::find(buffer.begin(), buffer.end(), state.start) !=
                buffer.end())
              size1 -= 1;
            buffer.clear();
            std::set_intersection(bagStore_[path_[state.index1 - 1]].begin(),
                                  bagStore_[path_[state.index1 - 1]].end(),
                                  bagStore_[path_[state.index2]].begin(),
                                  bagStore_[path_[state.index2]].end(),
                                  std::back_inserter(buffer));
            size2 = buffer.size();
            if (std::find(buffer.begin(), buffer.end(), state.start) !=
                buffer.end())
              size2 -= 1;
          }

          if (size2 > size1 || (width1 > width2 && size1 == size2)) {
            // move range2 before range1
            // not really an efficient way, but it works
            if (state.index2 != path_.size()) {
              std::vector<unsigned int> tmp(path_.begin() + state.index1,
                                            path_.end());
              unsigned int sep = state.index2 - state.index1;
              std::move(tmp.begin() + sep, tmp.end(),
                        path_.begin() + state.index1);
              std::move(tmp.begin(), tmp.begin() + sep,
                        path_.begin() + (path_.size() - sep));
            }
            state.index2 = state.index1 + (path_.size() - state.index2);
          }

          patches_.push_back(
              std::make_tuple(state.index1, state.index2 - 1, start));
          if (path_.size() != state.index2) {
            bagStore_[path_[state.index2]].insert(start);
            patches_.push_back(
                std::make_tuple(state.index2, path_.size() - 1, state.iend));
          }
        }
        newStart = state.iend;
      } else if (newStart == UINT32_MAX) {
        unsigned int v;
        if (std::distance(state.iter_start, state.iter_end) == 2) {
          unsigned int i = 0;
          if (state.iter_curr == state.iter_start) {
            state.index1 = path_.size();
            if (lastSubsumption != UINT32_MAX &&
                index[*state.iter_curr] != lastSubsumption)
              state.swapped = true;
          } else {
            state.index2 = path_.size();
            i = 1;
          }
          if (state.swapped)
            i = 1 - i;
          v = index[*(state.iter_start + i)];
        } else {
          if (state.iter_curr + 1 == state.iter_end) {
            addOnce = state.start;
          }
          v = index[*state.iter_curr];
        }
        state.iter_curr++;
        AddBag<1>({v});
        AddTask(v, state.iend);
        continue;
      }

      if (newStart != state.end) {
        if (handled_[newStart]) {
          stack_.pop_back();
          if (stack_.empty())
            return;
          continue;
        }
        handled_[newStart] = true;
      }
      start = newStart;

      auto iter = boost::adjacent_vertices(start, simplified_);
      state.iter_start = state.iter_curr = iter.first;
      state.iter_end = iter.second;
      state.iend = postdomtree_.GetIDom(start);
    }
  }

  void ApplyPatches(unsigned int index) {
    for (int p = patches_.size() - 1; p >= 0; p--) {
      auto &patch = patches_[p];
      unsigned int end = std::get<1>(patch);
      if (end < index)
        continue;
      unsigned int start = std::get<0>(patch);
      unsigned int v = std::get<2>(patch);

      for (unsigned int i = std::max(start, index); i < end; i++)
        bagStore_[path_[i]].insert(v);
      if (start < index) {
        std::get<1>(patch) = index - 1;
      } else {
        patches_.erase(patches_.begin() + p);
      }
    }
  }

  void AddTask(unsigned int start, unsigned int end) {
    if (start != end) {
      if (handled_[start])
        return;
      handled_[start] = true;
    }
    adjacency_iter_t iter_start, iter_end;
    boost::tie(iter_start, iter_end) =
        boost::adjacent_vertices(start, simplified_);
    stack_.push_back({iter_start, iter_start, iter_end, start, end,
                      postdomtree_.GetIDom(start), UINT32_MAX, UINT32_MAX,
                      false});
  }

  template <const int N> void AddBag(const std::array<unsigned int, N> &bag) {
    PathDecBag bag_;
    bag_.insert(bag.begin(), bag.end());
    if (addOnce != UINT32_MAX) {
      bag_.insert(addOnce);
      addOnce = UINT32_MAX;
    }
    if (path_.empty() || bagStore_[path_.back()] != bag_) {
      path_.push_back(bagStore_.size());
      bagStore_.push_back(bag_);
    }
  }

  // Try to turn the graph into a series-parallel graph.
  // Idea:
  // 1. Remove unreachable nodes.
  // 2. Remove loop backedges.
  // 3. Remove edges from loop body to outside.
  //
  // 1 is fine because we cannot reach those nodes from the entry node anyway.
  // For 2 and 3, we have to add these edges to a list of fixup, so we will
  // remember to add them back to the path decomposition later, to make sure
  // that the path decomposition is valid for the original graph. Also, after
  // removing an edge u->v, we have to check if u and v can still be reached
  // from the entry node, and can reach the exit node. For example, if u->v
  // corresponds to the instruction break, it is very likely that the out degree
  // of u is 0 after removing u->v, and can no longer reach the exit node. We
  // will add an edge from u to the exit node in this case. Another possible
  // case is an infinite loop with a break statement that corresponds to the
  // edge u->v, where v is the instruction after the loop. As the loop is an
  // infinite loop, the loop header will not connect to v, and removing u->v
  // will cause v to be unreachable from the entry node. In this case, we will
  // add the edge from loop header to v.
  // Also, it is possible for the loop to have no exit. In that case, we will
  // add an edge from the loop header to the exit node.
  //
  // Time complexity: O(|V|^2)
  // 1. Dominator tree calculation takes O(|V|^2) time.
  // 2. Backedge identification takes O(|V|^2) time. (Note: Probably not with
  // this implementation, we can use a hashset to make the Dominates check
  // constant time. But this implementation is faster in practice.)
  // 3. For each loop header, remove every edges from loop body to outside the
  // loop body. This takes O(|V|^2) time
  void RemoveBreakContinue() {
    auto index = boost::get(boost::vertex_index, simplified_);

    DomTree domtree(simplified_);

    lastVertex = boost::num_vertices(simplified_) - 1;
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
                          boost::container::small_vector<unsigned int, 20>>(
                    v, {u}));
          else
            (*entry).second.push_back(u);
        }
      }
    }

    boost::container::flat_map<unsigned int,
                               std::pair<PathDecBag, unsigned int>>
        loops;
    std::vector<unsigned int> loopSet(lastVertex + 1, UINT32_MAX);

    for (auto iter : backedges) {
      auto worklist =
          boost::container::small_vector<unsigned int, 20>(iter.second);
      // marked is the set of vertices inside the loop body
      auto marked = PathDecBag();
      marked.insert(iter.first);
      marked.insert(iter.second.begin(), iter.second.end());
      while (!worklist.empty()) {
        unsigned int u = worklist.back();
        worklist.pop_back();
        for (auto iter : boost::make_iterator_range(
                 boost::inv_adjacent_vertices(u, simplified_))) {
          unsigned int v = index[iter];
          if (marked.insert(v).second)
            worklist.push_back(v);
        }
      }
      marked.erase(iter.first);
      if (marked.empty())
        continue;
      unsigned int loopEnd = UINT32_MAX;
      // this is indeed a loop
      for (auto iter : boost::make_iterator_range(
               boost::adjacent_vertices(iter.first, simplified_))) {
        unsigned int v = index[iter];
        if (marked.find(v) == marked.end())
          loopEnd = v;
      }

      for (auto v : marked) {
        auto &l = loopSet[v];
        if (l == UINT32_MAX || loops[l].first.size() > marked.size())
          l = iter.first;
      }
      loops[iter.first].first = std::move(marked);
      loops[iter.first].second = loopEnd;
    }

    // the toRemove edges are left just to make postdomtree happy
    std::vector<std::pair<unsigned int, unsigned int>> toRemove;
    for (unsigned int u = 0; u <= lastVertex; u++) {
      auto loopHead = loopSet[u];
      if (loopHead == UINT32_MAX)
        continue;
      auto &marked = loops[loopHead];
      bool modified = true;
      while (modified) {
        modified = false;
        for (auto adj : boost::make_iterator_range(
                 boost::adjacent_vertices(u, simplified_))) {
          unsigned int v = index[adj];
          if (v == loopHead) {
            fixup_[u].insert(v);
            toRemove.push_back(std::make_pair(u, v));
          }
          if (v != marked.second && v != loopHead && marked.first.find(v) == marked.first.end()) {
            boost::remove_edge(u, v, simplified_);
            fixup_[u].insert(v);
            if (v != loopHead && !boost::edge(loopHead, v, simplified_).second) {
              boost::add_edge(loopHead, v, simplified_);
              if (marked.second == UINT32_MAX)
                marked.second = v;
            }
            if (boost::out_degree(u, simplified_) == 0) {
              unsigned int w = marked.second == UINT32_MAX ? loopHead : marked.second;
              boost::add_edge(u, w, simplified_);
              toRemove.push_back(std::make_pair(u, w));
            }
            // remove edge will invalidate iterator, so we have to break
            // the loop and loop again
            modified = true;
            break;
          }
        }
      }
    }
    for (auto &loop : loops) {
      if (loop.second.second == UINT32_MAX && loop.first != lastVertex) {
        // this is a loop with no exit
        // add an edge to the exit node so we can discover this loop
        boost::add_edge(loop.first, lastVertex, simplified_);
        toRemove.push_back(std::make_pair(loop.first, lastVertex));
      }
    }
    auto tmp = boost::make_reverse_graph(simplified_);
    postdomtree_ = DomTree(tmp, lastVertex);

    for (auto iter : toRemove) {
      boost::remove_edge(iter.first, iter.second, simplified_);
    }
  }
};

