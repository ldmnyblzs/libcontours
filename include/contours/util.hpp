#ifndef UTIL_HPP
#define UTIL_HPP 1

#include <CGAL/Kernel_traits.h>
#include <CGAL/boost/graph/iterator.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/range/join.hpp>

namespace contours {
  template <typename Graph>
  auto every_adjacent(const typename boost::graph_traits<Graph>::vertex_descriptor vertex, const Graph &graph) {
    return boost::join(boost::adjacent_vertices(vertex, graph), boost::inv_adjacent_vertices(vertex, graph));
  }
template <typename Mesh, typename PointMap>
auto edge_segment(
    const Mesh &mesh, const PointMap &point,
    const typename boost::graph_traits<Mesh>::edge_descriptor edge) {
  using namespace boost;
  using namespace CGAL;
  using Point = typename property_traits<PointMap>::value_type;
  using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
  using Segment = typename Kernel::Segment_3;

  const auto h = halfedge(edge, mesh);
  return Segment(get(point, source(h, mesh)), get(point, target(h, mesh)));
}

template <typename Mesh, typename PointMap>
auto face_triangle(
    const Mesh &mesh, const PointMap &point,
    const typename boost::graph_traits<Mesh>::face_descriptor face) {
  using namespace boost;
  using namespace CGAL;
  using Point = typename property_traits<PointMap>::value_type;
  using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
  using Triangle = typename Kernel::Triangle_3;

  const auto h = halfedge(face, mesh);
  return Triangle(get(point, source(h, mesh)), get(point, target(h, mesh)),
                  get(point, target(next(h, mesh), mesh)));
}

  template <typename Vector>
  inline auto normalize(const Vector &vector) {
    return vector / sqrt(vector.squared_length());
  }

  template <typename Graph>
  auto in_edge(const typename boost::graph_traits<Graph>::vertex_descriptor vertex,
	       const Graph &graph) {
    using namespace boost;
    return *begin(in_edges(vertex, graph));
  }
  template <typename Graph>
  auto out_edge(const typename boost::graph_traits<Graph>::vertex_descriptor vertex,
	       const Graph &graph) {
    using namespace boost;
    return *begin(out_edges(vertex, graph));
  }
  template <typename Graph>
  auto in_source(const typename boost::graph_traits<Graph>::vertex_descriptor vertex,
		 const Graph &graph) {
    return boost::source(in_edge(vertex, graph), graph);
  }
  template <typename Graph>
  auto out_target(const typename boost::graph_traits<Graph>::vertex_descriptor vertex,
		 const Graph &graph) {
    return boost::target(out_edge(vertex, graph), graph);
  }
}; // namespace shape

#endif // UTIL_HPP
