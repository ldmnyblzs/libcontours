#ifndef MESH_PROPERTIES_HPP
#define MESH_PROPERTIES_HPP

#include <CGAL/Kernel/global_functions_3.h>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/algorithm/min_element.hpp>
#include <boost/range/numeric.hpp>
#include <cmath>
#include <functional>

#include "util.hpp"

namespace contours {
template <typename Mesh, typename PointMap, typename Point>
auto max_distance(const Mesh &mesh, const PointMap &point,
                  const Point &origin) {
  using namespace std;
  using namespace boost;
  using namespace boost::adaptors;
  using namespace CGAL;
  using std::ref;

  const auto transformer = [&point, &origin](const auto &v) {
    return squared_distance(origin, get(point, v));
  };
  return sqrt(*max_element(vertices(mesh) | transformed(ref(transformer))));
}

template <typename Mesh, typename PointMap, typename Point>
auto min_distance(const Mesh &mesh, const PointMap &point,
                  const Point &origin) {
  using namespace std;
  using namespace boost;
  using namespace boost::adaptors;
  using namespace CGAL;
  using std::ref;

  const auto transformer = [&mesh, &point, &origin](const auto &f) {
    return squared_distance(origin, face_triangle(mesh, point, f));
  };
  return sqrt(*min_element(faces(mesh) | transformed(ref(transformer))));
}

template <typename Mesh, typename PointMap>
auto surface_area(const Mesh &mesh, const PointMap &point) {
  using namespace boost;
  return accumulate(
      faces(mesh), 0.0, [&mesh, &point](const auto init, const auto face) {
        return init + sqrt(face_triangle(mesh, point, face).squared_area());
      });
}

template <typename Mesh, typename PointMap, typename Point>
auto volume(const Mesh &mesh, const PointMap &point, const Point &centroid) {
  using namespace boost;
  using namespace CGAL;
  return accumulate(faces(mesh), 0.0, [&](const auto init, const auto face) {
    const auto h = halfedge(face, mesh);
    return init + abs(volume(centroid, get(point, source(h, mesh)),
                             get(point, target(h, mesh)),
                             get(point, target(next(h, mesh), mesh))));
  });
}
}; // namespace shape

#endif // MESH_PROPERTIES_HPP
