#ifndef MAKE_REEB_HPP
#define MAKE_REEB_HPP 1

#include <boost/container/flat_set.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <iostream>
#include <vector>

namespace shape {
template <typename Graph, typename Edges, typename RemoveMap>
static void
mark_to(Graph &graph,
        const typename boost::graph_traits<Graph>::vertex_descriptor root,
        const Edges &edges, RemoveMap &remove) {
  using namespace boost;
  auto degree = 0;
  for (const auto adjacent : make_iterator_range(in_edges(root, graph)))
    if (!get(remove, source(adjacent, graph)))
      degree++;
  if (degree <= 1) {
    for (const auto edge : make_iterator_range(out_edges(root, graph)))
      if (find(edges, edge) == edges.end())
        mark_to(graph, target(edge, graph), edges, remove);
      else
        return;
    put(remove, root, true);
  }
}
template <typename Graph, typename Roots, typename Edges, typename RemoveMap>
void mark_inside(Graph &graph, const Roots &roots, const Edges &edges,
                 RemoveMap &remove) {
  for (const auto r : roots)
    mark_to(graph, r, edges, remove);
}

template <typename Graph, typename RemoveMap, typename EdgeLevelMap,
          typename VertexLevelMap>
void make_reeb(Graph &graph, const RemoveMap &remove,
               const EdgeLevelMap &edge_level, VertexLevelMap &vertex_level) {
  using namespace std;
  using namespace boost;
  using namespace boost::adaptors;
  using Vertex = typename graph_traits<Graph>::vertex_descriptor;

  // remove marked vertices
  vector<Vertex> to_remove;
  for (const auto vertex : make_iterator_range(vertices(graph))) {
    if (get(remove, vertex)) {
      clear_vertex(vertex, graph);
      to_remove.push_back(vertex);
    }
  }
  for (const auto vertex : reverse(to_remove))
    remove_vertex(vertex, graph);

  // set vertex level for equilibria
  for (const auto vertex : make_iterator_range(vertices(graph))) {
    const auto din = in_degree(vertex, graph);
    const auto dout = out_degree(vertex, graph);

    if (din == 0) {
      put(vertex_level, vertex, get(edge_level, out_edge(vertex, graph)) - 1);
    } else if (dout == 0) {
      put(vertex_level, vertex, get(edge_level, in_edge(vertex, graph)) + 1);
    } else if (din > 1 || dout > 1) {
      put(vertex_level, vertex, get(edge_level, in_edge(vertex, graph)) + 0.5);
    }
  }

  // isolate and remove vertices of degree 2
  to_remove.clear();
  for (const auto vertex : make_iterator_range(vertices(graph))) {
    const auto din = in_degree(vertex, graph);
    const auto dout = out_degree(vertex, graph);

    if (din == 1 && dout == 1) {
      add_edge(source(in_edge(vertex, graph), graph),
               target(out_edge(vertex, graph), graph), graph);
      clear_vertex(vertex, graph);
      to_remove.push_back(vertex);
    }
  }
  for (const auto vertex : reverse(to_remove))
    remove_vertex(vertex, graph);
}

template <typename Graph, typename VertexLevelMap, typename VertexIndex>
bool make_morse(const Graph &graph, const VertexLevelMap &level,
                VertexIndex &index) {
  using boost::container::flat_set;
  using Real = typename boost::property_traits<VertexLevelMap>::value_type;
  using namespace boost;

  flat_set<Real> levels;
  for (const auto vertex : make_iterator_range(vertices(graph)))
    if (!levels.insert(get(level, vertex)).second)
      return false;
  for (const auto vertex : make_iterator_range(vertices(graph)))
    put(index, vertex, levels.index_of(levels.find(get(level, vertex))));
  return true;
}
}; // namespace shape

#endif // MAKE_REEB_HPP
