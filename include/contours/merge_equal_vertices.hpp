/*
  Copyright 2019 Balázs Ludmány

  This file is part of libcontours.

  libcontours is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  libcontours is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with libcontours.  If not, see <https://www.gnu.org/licenses/>.
*/

/*
  Copyright 2019 Balázs Ludmány

  This file is part of libcontours.

  libcontours is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  libcontours is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with libcontours.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef MERGE_EQUAL_VERTICES_HPP
#define MERGE_EQUAL_VERTICES_HPP 1

#include <boost/graph/graph_traits.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/iterator_range.hpp>
#include <deque>
#include <vector>
#include <boost/range/algorithm/find.hpp>
#include <boost/property_map/property_map.hpp>

namespace contours {
  template <typename Graph, typename ArcListMap, typename EdgeProperty>
static void
add_edge_unique(const typename boost::graph_traits<Graph>::vertex_descriptor &source,
                const typename boost::graph_traits<Graph>::vertex_descriptor &target,
		ArcListMap &arc_list,
		typename boost::property_traits<ArcListMap>::value_type &other_arcs,
		const EdgeProperty &property,
                Graph &graph) {
  const auto existing = edge(source, target, graph);
  if (existing.second)
    arc_list[existing.first].splice(arc_list[existing.first].cend(), other_arcs);
  else
    add_edge(source, target, property, graph);
}

template <typename Graph, typename EqualVerticesMap, typename AreaMap,
          typename VisitedMap, typename ArcListMap>
void merge_equal_vertices(Graph &graph, EqualVerticesMap &equal, AreaMap &area,
                          VisitedMap &visited, ArcListMap &arc_list) {
  using namespace std;
  using namespace boost;
  using namespace boost::adaptors;
  using Vertex = typename graph_traits<Graph>::vertex_descriptor;
  
  for (const auto vertex : make_iterator_range(vertices(graph))) {
    if (!get(visited, vertex)) {
      visited[vertex] = true;
      deque<Vertex> equal_vertices(equal[vertex].cbegin(), equal[vertex].cend());

      while (!equal_vertices.empty()) {
        if (!get(visited, equal_vertices.front())) {
          visited[equal_vertices.front()] = true;

          for (const auto edge :
               make_iterator_range(out_edges(equal_vertices.front(), graph)))
            add_edge_unique(vertex, target(edge, graph), arc_list, get(arc_list, edge), get(edge_all, graph, edge), graph);
          for (const auto edge :
               make_iterator_range(in_edges(equal_vertices.front(), graph)))
            add_edge_unique(source(edge, graph), vertex, arc_list, get(arc_list, edge), get(edge_all, graph, edge), graph);

	  clear_vertex(equal_vertices.front(), graph);
          area[vertex] += get(area, equal_vertices.front());
          for (const auto eq : equal[equal_vertices.front()])
	    if (!get(visited, eq))
	      equal_vertices.push_back(eq);
        }
	equal_vertices.pop_front();
      }
    }
  }

  vector<Vertex> to_remove;
  for (const auto vertex : make_iterator_range(vertices(graph))) {
    if (degree(vertex, graph) == 0)
      to_remove.push_back(vertex);
  }

  for (const auto vertex : reverse(to_remove))
    remove_vertex(vertex, graph);
}
}; // namespace shape

#endif // MERGE_EQUAL_VERTICES_HPP
