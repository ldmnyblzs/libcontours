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

#ifndef INTERSECT_EDGES_HPP
#define INTERSECT_EDGES_HPP

#include <CGAL/Origin.h>
#include <boost/property_map/property_map.hpp>
#include <boost/range/irange.hpp>
#include <map>

#include "util.hpp"

namespace contours {
template <typename Mesh, typename PointMap, typename Point, typename Real,
          typename HalfedgeIntersectionMap>
void intersect_halfedges(const Mesh &mesh, const PointMap &point,
                         const Point &origin, const Real min_distance,
                         const Real step, HalfedgeIntersectionMap &intersection) {
  using namespace std;
  using namespace CGAL;
  using namespace boost;
  using CGAL::max;
  using Intersection =
      typename property_traits<HalfedgeIntersectionMap>::value_type::value_type;

  for (const auto edge : edges(mesh)) {
    const auto segment = edge_segment(mesh, point, edge);
    const auto edge_vector = segment.to_vector();
    const auto min_edge_distance = sqrt(squared_distance(origin, segment));
    const auto source_sq_distance = squared_distance(origin, segment[0]);
    const auto max_edge_distance =
        sqrt(max(source_sq_distance, squared_distance(origin, segment[1])));

    const int first_level =
        floor((min_edge_distance - min_distance) / step) + 1;
    const int next_level = ceil((max_edge_distance - min_distance) / step);

    if (first_level < next_level) {
      const auto h = halfedge(edge, mesh);
      auto &map = get(intersection, h);
      auto &opposite_map = get(intersection, opposite(h, mesh));

      const auto a = segment.squared_length();
      const auto b = 2.0 * scalar_product(edge_vector, segment[0] - origin);

      for (const int level : irange(first_level, next_level)) {
        const Real r = min_distance + level * step;
        const auto c = source_sq_distance - r * r;
        const auto D = b * b - 4.0 * a * c;
	const auto u1 = (-b - sqrt(D)) / (2.0 * a);
	const bool has_u1 = u1 > 0.0 && u1 < 1.0;
	const auto u2 = (-b + sqrt(D)) / (2.0 * a);
	const bool has_u2 = u2 > 0.0 && u2 < 1.0;
	if (has_u1 && has_u2) {
	  const auto I1 = segment[0] + u1 * edge_vector;
	  const auto I2 = segment[0] + u2 * edge_vector;
	  map.emplace_hint(map.end(), piecewise_construct,
			   forward_as_tuple(level), forward_as_tuple(I1, I2));
	  opposite_map.emplace_hint(opposite_map.end(), piecewise_construct,
				    forward_as_tuple(level),
				    forward_as_tuple(I2, I1));
	} else if (has_u1) {
	  const auto I1 = segment[0] + u1 * edge_vector;
	  const Intersection p(piecewise_construct,
			       forward_as_tuple(level),
			       forward_as_tuple(I1, I1));
	  map.insert(map.end(), p);
	  opposite_map.insert(opposite_map.end(), p);
	} else if (has_u2) {
	  const auto I2 = segment[0] + u2 * edge_vector;
	  const Intersection p(piecewise_construct,
			       forward_as_tuple(level),
			       forward_as_tuple(I2, I2));
	  map.insert(map.end(), p);
	  opposite_map.insert(opposite_map.end(), p);
	}
      }
    }
  }
}
}; // namespace shape

#endif // INTERSECT_EDGES_HPP
