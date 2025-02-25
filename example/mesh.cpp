#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/centroid.h>
#include <boost/container/flat_set.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <iostream>

#include "contours/discover_graph.hpp"
#include "contours/find_equilibra.hpp"
#include "contours/intersect_faces.hpp"
#include "contours/intersect_halfedges.hpp"
#include "contours/merge_equal_vertices.hpp"
#include "contours/mesh_properties.hpp"
#include "contours/stl_io.hpp"
#include "contours/make_reeb.hpp"
#include "contours/encode_graph.hpp"

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Vector = Kernel::Vector_3;
using Mesh = CGAL::Surface_mesh<Point>;
using Halfedge = Mesh::Halfedge_index;

struct VertexProperty;
struct EdgeProperty;

using Graph =
    boost::adjacency_list<boost::vecS, boost::listS, boost::bidirectionalS,
                          VertexProperty, EdgeProperty>;
using Reverse = boost::reverse_graph<Graph>;

using GraphVertex = boost::graph_traits<Graph>::vertex_descriptor;
using GraphEdge = boost::graph_traits<Graph>::edge_descriptor;
using ReverseVertex = boost::graph_traits<Reverse>::vertex_descriptor;
using ReverseEdge = boost::graph_traits<Reverse>::edge_descriptor;

struct VertexProperty {
  double area = 0.0;
  boost::container::flat_set<GraphVertex> eq_edges;
  bool visited = false;
  int id;
  float level;
  int label;
};

struct Arc {
  Point center;
  Point source;
  Point target;
  Vector normal;
  Arc(Point center, Point source, Point target, Vector normal) :
    center(std::move(center)),
    source(std::move(source)),
    target(std::move(target)),
    normal(std::move(normal)) {
  }
};

struct EdgeProperty {
  double area_inside = 0.0, area_outside = 0.0;
  boost::container::flat_set<GraphVertex> roots_inside, roots_outside;
  int level = 0;
  std::list<Arc> arcs;
};

int main(int argc, char **argv) {
  Mesh mesh;
  Graph graph;

  const boost::iostreams::mapped_file_source file(argv[1]);
  contours::read_STL(file.begin(), file.end(), mesh, mesh.points());
  const int count = atoi(argv[2]);
  const double ratio = atof(argv[3]);

  const Point origin =
      CGAL::centroid(mesh.points().begin(), mesh.points().end());
  const double area = contours::surface_area(mesh, mesh.points());
  const double max_distance = contours::max_distance(mesh, mesh.points(), origin);
  const double min_distance = contours::min_distance(mesh, mesh.points(), origin);
  const double step = (max_distance - min_distance) / (count + 1);

  auto halfedge_intersections =
      mesh.add_property_map<Halfedge, std::map<int, std::pair<Point, Point>>>(
              "h:intersections")
          .first;

  contours::intersect_halfedges(mesh, mesh.points(), origin, min_distance, step,
                             halfedge_intersections);

  auto to_halfedge =
      mesh.add_property_map<Halfedge, std::map<int, GraphVertex>>(
              "h:to_halfedge")
          .first;
  auto from_halfedge =
      mesh.add_property_map<Halfedge, std::map<int, GraphVertex>>(
              "h:from_halfedge")
          .first;
  auto area_map = boost::get(&VertexProperty::area, graph);
  auto eq_map = boost::get(&VertexProperty::eq_edges, graph);
  auto edge_level = boost::get(&EdgeProperty::level, graph);
  auto arc_list = boost::get(&EdgeProperty::arcs, graph);
  contours::intersect_faces(mesh, mesh.points(), halfedge_intersections,
                         to_halfedge, from_halfedge, origin, min_distance, step,
                         graph, area_map, eq_map, edge_level, arc_list);

  auto visited_map = boost::get(&VertexProperty::visited, graph);
  contours::merge_equal_vertices(graph, eq_map, area_map, visited_map, arc_list);

  auto area_inside_map = boost::get(&EdgeProperty::area_inside, graph);
  auto roots_inside_map = boost::get(&EdgeProperty::roots_inside, graph);
  std::vector<GraphVertex> stable_vertices, unstable_vertices;
  contours::discover_graph(graph, area_map, area_inside_map, roots_inside_map,
                        std::back_inserter(stable_vertices));

  auto reverse = make_reverse_graph(graph);
  auto area_outside_map = boost::get(&EdgeProperty::area_outside, reverse);
  auto roots_outside_map = boost::get(&EdgeProperty::roots_outside, reverse);
  contours::discover_graph(reverse, area_map, area_outside_map, roots_outside_map,
                        std::back_inserter(unstable_vertices));

  std::vector<GraphEdge> stable_edges;
  std::vector<ReverseEdge> unstable_edges;
  contours::find_equilibria(graph, stable_vertices, area_inside_map,
                         roots_inside_map, std::back_inserter(stable_edges),
                         area * ratio);
  contours::find_equilibria(reverse, unstable_vertices, area_outside_map,
                         roots_outside_map,
                         std::back_inserter(unstable_edges),
                         area * ratio);
  std::cout << "S=" << stable_edges.size() << " U=" << unstable_edges.size() << '\n';

  for (const auto vertex : boost::make_iterator_range(boost::vertices(graph)))
    graph[vertex].visited = false;
  
  contours::mark_inside(graph, stable_vertices, stable_edges, visited_map);
  contours::mark_inside(reverse, unstable_vertices, unstable_edges, visited_map);

  edge_level = boost::get(&EdgeProperty::level, graph);
  auto vertex_level = boost::get(&VertexProperty::level, graph);
  contours::make_reeb(graph, visited_map, edge_level, vertex_level);
  auto vertex_index = boost::get(&VertexProperty::id, graph);
  int ind = 0;
  for (const auto vertex : boost::make_iterator_range(boost::vertices(graph)))
    graph[vertex].id = ind++;
  auto vertex_label = boost::get(&VertexProperty::label, graph);
  contours::reeb_encode(graph, vertex_index, vertex_label);
  std::cout << "Reeb: " << contours::encode(graph, vertex_label) << '\n';

  std::ofstream graph_file("out.dot");
  boost::write_graphviz(
      graph_file, graph,
      boost::make_label_writer(vertex_label),
      boost::default_writer(),
      boost::default_writer(),
      get(&VertexProperty::id, graph));
  
  if (contours::make_morse(graph, vertex_level, vertex_label))
    std::cout << "Morse: " << contours::encode(graph, vertex_label) << '\n';
}
