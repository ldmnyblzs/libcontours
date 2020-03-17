#ifndef ENCODE_GRAPH_HPP
#define ENCODE_GRAPH_HPP 1

#include <bliss-0.73/graph.hh>
#include <boost/container/static_vector.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/multi_array.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/range/algorithm/generate.hpp>
#include <boost/range/algorithm/permutation.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/combine.hpp>
#include <deque>
#include <map>
#include <sstream>
#include <tuple>
#include <vector>

#include "util.hpp"

namespace contours {
template <typename Graph, typename IndexMap>
std::string encode(const Graph &graph, const IndexMap &index) {
  using namespace std;
  using namespace boost;
  // using namespace boost::algorithm;
  using Id = typename boost::property_traits<IndexMap>::value_type;

  vector<pair<Id, Id>> endpoints;
  for (const auto edge : make_iterator_range(edges(graph))) {
    endpoints.emplace_back(get(index, source(edge, graph)),
                           get(index, target(edge, graph)));
  }
  sort(endpoints);
  stringstream encoded;
  for (const auto e : endpoints) {
    encoded << hex << e.first << e.second;
  }
  return encoded.str();
}

template <typename Graph, typename IndexMap>
bool is_valid(const Graph &graph, const IndexMap &index) {
  using namespace boost;

  int i = 0;
  for (const auto vertex : make_iterator_range(vertices(graph))) {
    const auto vertex_index = get(index, vertex);
    for (const auto in_vertex :
         make_iterator_range(inv_adjacent_vertices(vertex, graph)))
      if (get(index, in_vertex) > vertex_index)
        return false;
    for (const auto out_vertex :
         make_iterator_range(adjacent_vertices(vertex, graph)))
      if (get(index, out_vertex) < vertex_index)
        return false;
  }
  return true;
}
template <typename Graph, typename IndexMap, typename LabelMap>
void reeb_encode(const Graph &graph, const IndexMap &index, LabelMap &label) {
  using namespace bliss;
  using namespace boost;

  Digraph digraph(num_vertices(graph));
  for (const auto edge : make_iterator_range(edges(graph)))
    digraph.add_edge(get(index, source(edge, graph)),
                      get(index, target(edge, graph)));

  Stats stats;
  const unsigned int *canonical =
      digraph.canonical_form(stats, nullptr, nullptr);

  for (const auto vertex : make_iterator_range(vertices(graph)))
    put(label, vertex, canonical[get(index, vertex)]);
}
}; // namespace shape

#endif // ENCODE_GRAPH_HPP
