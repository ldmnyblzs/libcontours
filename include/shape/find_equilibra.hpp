#ifndef FIND_EQUILIBRIA_HPP
#define FIND_EQUILIBRIA_HPP 1

#include <boost/container/flat_set.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/remove_if.hpp>
#include <boost/range/algorithm/set_algorithm.hpp>
#include <optional>
#include <utility>

#include "util.hpp"

namespace shape {
template <typename Graph, typename AreaMap, typename Real>
static std::optional<typename boost::graph_traits<Graph>::edge_descriptor>
search_bigger(const Graph &graph,
              const typename boost::graph_traits<Graph>::vertex_descriptor root,
              const AreaMap &edge_area, const Real area) {
  using namespace std;
  if (out_degree(root, graph) != 1)
    return nullopt;

  const auto next = *out_edges(root, graph).first;
  if (get(edge_area, next) >= area)
    return next;

  return search_bigger(graph, target(next, graph), edge_area, area);
}

template <typename Graph, typename VertexRange, typename AreaMap,
          typename RootsMap, typename EdgeIterator, typename Real>
void find_equilibria(const Graph &graph, const VertexRange &vertices,
                     AreaMap &edge_area, RootsMap &roots, EdgeIterator edges,
                     Real area) {
  using namespace boost;
  using boost::container::flat_set;
  using namespace boost::range;
  using std::ref;

  using Edge = typename graph_traits<Graph>::edge_descriptor;

  flat_set<Edge> edge_set;

  for (const auto &vertex : vertices)
    if (const auto edge = search_bigger(graph, vertex, edge_area, area))
      edge_set.insert(*edge);

  edge_set.erase(remove_if(edge_set,
                           [&edge_set, &roots](const auto edge) noexcept {
                             for (const auto other : edge_set)
                               if (edge != other && includes(get(roots, edge),
                                                             get(roots, other)))
                                 return true;
                             return false;
                           }),
                 edge_set.end());

  copy(edge_set, edges);
}
}; // namespace shape

#endif // FIND_EQUILIBRIA_HPP
