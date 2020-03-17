#ifndef AXES_HPP
#define AXES_HPP

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/squared_distance_3.h>
#include <array>
#include <boost/property_map/property_map.hpp>

#include "util.hpp"

namespace contours {

template <typename PointMap>
std::array<typename CGAL::Kernel_traits<
	     typename boost::property_traits<PointMap>::value_type>::Kernel::Vector_3,
           3>
axes(const PointMap &points) {
  using Point = typename boost::property_traits<PointMap>::value_type;
  using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
  using Vector = typename Kernel::Vector_3;
  using Plane = typename Kernel::Plane_3;
  using Line = typename Kernel::Line_3;

  double max_squared_a = 0, max_squared_b = 0, max_squared_c = 0;
  Vector a, b, c;

  for (auto point1 = points.begin(); point1 != std::prev(points.end());
       ++point1) {
    for (auto point2 = std::next(point1); point2 != points.end(); ++point2) {
      const auto squared_distance = CGAL::squared_distance(*point1, *point2);
      if (squared_distance > max_squared_a) {
        max_squared_a = squared_distance;
        a = Vector(*point1, *point2);
      }
    }
  }

  const Plane projection_plane(CGAL::ORIGIN, a);
  for (auto point1 = points.begin(); point1 != std::prev(points.end());
       ++point1) {
    const auto projected1 = projection_plane.projection(*point1);
    for (auto point2 = std::next(point1); point2 != points.end(); ++point2) {
      const auto projected2 = projection_plane.projection(*point2);
      const auto squared_distance =
          CGAL::squared_distance(projected1, projected2);
      if (squared_distance > max_squared_b) {
        max_squared_b = squared_distance;
        b = Vector(projected1, projected2);
      }
    }
  }

  c = CGAL::cross_product(a, b);
  const Line projection_line(CGAL::ORIGIN, c);
  for (auto point1 = points.begin(); point1 != std::prev(points.end());
       ++point1) {
    const auto projected1 = projection_line.projection(*point1);
    for (auto point2 = std::next(point1); point2 != points.end(); ++point2) {
      const auto projected2 = projection_line.projection(*point2);
      const auto squared_distance =
          CGAL::squared_distance(projected1, projected2);
      if (squared_distance > max_squared_c)
        max_squared_c = squared_distance;
    }
  }

  c /= CGAL::sqrt(c.squared_length());
  c *= CGAL::sqrt(max_squared_c);

  return {a, b, c};
}

template <typename PointMap, typename Vector>
std::pair<double, double> projected_properties(const PointMap &points,
                                               const Vector &a, const Vector &b,
                                               const Vector &c) {
  using Point = typename boost::property_traits<PointMap>::value_type;
  using Kernel = typename CGAL::Kernel_traits<Vector>::Kernel;
  using Plane = typename Kernel::Plane_3;

  const auto x = normalize(a);
  const auto y = normalize(b);
  const auto z = normalize(c);

  const Plane plane(CGAL::ORIGIN, c);
  const CGAL::Aff_transformation_3<Kernel> new_basis(
      x.x(), x.y(), x.z(), y.x(), y.y(), y.z(), z.x(), z.y(), z.z());

  std::vector<typename Kernel::Point_2> projected_points;
  projected_points.reserve(boost::size(points));

  for (const auto &point : points) {
    const auto projected = new_basis(plane.projection(point));
    projected_points.emplace_back(projected.x(), projected.y());
  }

  CGAL::Polygon_2<Kernel> outline;
  CGAL::convex_hull_2(projected_points.cbegin(), projected_points.cend(),
                      std::back_inserter(outline));

  const auto circumference =
      std::accumulate(outline.edges_begin(), outline.edges_end(), 0.0,
                      [](const auto init, const auto &edge) {
                        return init + CGAL::sqrt(edge.squared_length());
                      });
  return {circumference, outline.area()};
}
}; // namespace shape

#endif // AXES_HPP
