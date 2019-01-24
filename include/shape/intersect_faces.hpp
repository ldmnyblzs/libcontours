#ifndef INTERSECT_FACES_HPP
#define INTERSECT_FACES_HPP

#include <CGAL/Kernel_traits.h>
#include <bitset>
#include <boost/container/static_vector.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/irange.hpp>
#include <map>
#include <optional>
#include "util.hpp"

namespace shape {
template <typename Point, typename Vector>
static auto arc_area(const Point &source, const Point &target,
                     const Vector &normal, const Point &center) {
  const auto to_source = source - center;
  const auto to_target = target - center;
  const auto angle =
      atan2(-scalar_product(normal, cross_product(to_source, to_target)),
            -scalar_product(to_source, to_target)) +
      M_PI;
  return 0.5 * to_source.squared_length() * (angle - sin(angle));
}

template <typename Point, typename Vector>
static auto polygon_area(const std::vector<Point> &points,
                         const Vector &normal, const Point &origin) {
  using namespace CGAL;
  const auto size = points.size();
  if (size >= 3) {
    Vector sum(0, 0, 0);
    for (std::size_t j = 0; j < size; ++j)
      sum += cross_product(points[j] - origin, points[(j + 1) % size] - origin);
    return scalar_product(sum, normal) * 0.5;
  }
  return 0.0;
}

static bool is_subset(const std::bitset<3> &superset,
                      const std::bitset<3> &subset) {
  return (superset & subset) == subset;
}

template <typename GraphVertex, typename EquivalentMap, typename Graph>
static void add_edge_between_faces(std::map<int, GraphVertex> &lookup_map,
                                   std::map<int, GraphVertex> &insert_map,
                                   int level, const GraphVertex vertex,
                                   EquivalentMap &equivalent, Graph &graph) {
  const auto it = lookup_map.find(level);
  if (it != lookup_map.end()) {
    equivalent[vertex].insert(it->second);
    equivalent[it->second].insert(vertex);
    lookup_map.erase(it);
  } else {
    insert_map.emplace(level, vertex);
  }
}

static constexpr std::array<std::bitset<3>, 3> edge_masks = {0b101, 0b011,
                                                             0b110};

template <typename Mesh, typename PointMap, typename EdgeIntersectionMap,
          typename HalfedgeVertexMap, typename Point, typename Real,
          typename Graph, typename VertexAreaMap, typename EquivalentMap,
          typename LevelMap, typename ArcListMap>
void intersect_faces(const Mesh &mesh, const PointMap &point,
                     const EdgeIntersectionMap &intersections,
                     HalfedgeVertexMap &to_halfedge,
                     HalfedgeVertexMap &from_halfedge, const Point &origin,
                     const Real min_distance, const Real step, Graph &graph,
                     VertexAreaMap &area, EquivalentMap &equivalent,
                     LevelMap &edge_level, ArcListMap &arc_list) {
  using std::array;
  using std::move;
  using std::optional;
  using namespace std;
  using namespace boost;
  using container::static_vector;
  using namespace CGAL;
  using CGAL::max;
  using CGAL::min;
  using Kernel = typename Kernel_traits<Point>::Kernel;
  using Sphere = typename Kernel::Sphere_3;
  using Circle = typename Kernel::Circle_3;
  using Halfedge = typename graph_traits<Mesh>::halfedge_descriptor;
  using iterator =
      typename property_traits<EdgeIntersectionMap>::value_type::const_iterator;
  using GraphVertex = typename graph_traits<Graph>::vertex_descriptor;
  using ArcList = typename property_traits<ArcListMap>::value_type;
  using Arc = typename ArcList::value_type;

  for (const auto face : faces(mesh)) {
    // get random access to a few things around the face
    const auto triangle = face_triangle(mesh, point, face);
    const auto face_plane = triangle.supporting_plane();
    const auto face_normal = normalize(face_plane.orthogonal_vector());
    const auto face_center = face_plane.projection(origin);
    static_vector<Halfedge, 3> halfedges, opposites;
    static_vector<iterator, 3> begins, ends;
    static_vector<double, 3> source_distances;
    int min_level = numeric_limits<int>::max();
    int next_level = numeric_limits<int>::min();
    int empty_edges = 0;
    for (const auto h : halfedges_around_face(halfedge(face, mesh), mesh)) {
      halfedges.push_back(h);
      opposites.push_back(opposite(h, mesh));
      const auto &inter = intersections[h];
      begins.push_back(inter.cbegin());
      ends.push_back(inter.cend());
      if (inter.empty()) {
        ++empty_edges;
      } else {
        next_level = max(next_level, prev(ends.back())->first + 1);
        min_level = min(min_level, begins.back()->first);
      }
      source_distances.push_back(
          sqrt(squared_distance(origin, triangle[halfedges.size() - 1])));
    }

    const int first_level =
        floor((sqrt(squared_distance(origin, triangle)) - min_distance) /
              step) +
        1;
    if (empty_edges == 3) {
      const auto max_distance = max_element(source_distances);
      next_level = min_level = ceil((*max_distance - min_distance) / step);
    }

    vector<pair<Point, Point>> previous_arcs;
    vector<bitset<3>> previous_covers;
    vector<pair<int, int>> previous_edges;
    vector<GraphVertex> previous_vertices;

    previous_vertices.push_back(add_vertex(graph));

    // [first_level:min_level): circles
    for (const int level : irange(first_level, min_level)) {
      const Real r = min_distance + level * step;
      const auto new_vertex = add_vertex(graph);
      const auto new_edge = add_edge(previous_vertices[0], new_vertex, graph).first;
      put(edge_level, new_edge, level);

      Circle circle(face_center, r * r - squared_distance(origin, face_center), face_normal);
      const double circle_area = circle.approximate_area();
      area[previous_vertices[0]] += circle_area;
      area[new_vertex] -= circle_area;
      previous_vertices[0] = new_vertex;
      const auto base = normalize(face_plane.base1()) * sqrt(circle.squared_radius());
      const auto p1 = face_center + base;
      const auto p2 = face_center - base;
      const ArcList arcs{Arc(face_center, p1, p2, face_normal), Arc(face_center, p2, p1, face_normal)};
      put(arc_list, new_edge, arcs);
    }

    previous_covers.push_back(0b111);

    // [min_level:next_level): arcs
    for (const int level : irange(min_level, next_level)) {
      const Real r = min_distance + level * step;

      vector<pair<Point, Point>> current_arcs;
      vector<bitset<3>> current_covers;
      vector<pair<int, int>> current_edges;
      vector<GraphVertex> current_vertices;

      bitset<3> passed_by;
      int last_edge = 2;
      optional<Point> previous_point;
      int previous_edge;

      // find current arcs
      for (int e = 0; e <= last_edge; ++e) {
        const auto edge = e % 3;
        if (begins[edge] != ends[edge] && begins[edge]->first == level) {
          if (previous_point && source_distances[edge] > r) {
            current_arcs.emplace_back(*previous_point,
                                      begins[edge]->second.first);
            current_covers.push_back(passed_by);
            current_edges.emplace_back(previous_edge, edge);
            current_vertices.push_back(add_vertex(graph));
          }

          passed_by.reset();
          previous_point = begins[edge]->second.second;
          previous_edge = edge;

          if (last_edge == 2)
            last_edge = 3 + edge;
          else
            ++begins[edge];
        }
        passed_by.set(edge);
      }

      // match up current arcs with previous_arcs
      for (const auto p : irange(0ul, previous_covers.size())) {
        const auto to_cover = previous_covers.at(p).count();
        array<int, 3> covers{3, 3, 3};

        for (const auto c : irange(0ul, current_covers.size())) {
          if (is_subset(previous_covers.at(p), current_covers.at(c))) {
            previous_covers.at(p) &= ~current_covers.at(c);

            bool first = true;
            for (const int i : irange(0, 3)) {
              if (current_covers.at(c)[i]) {
                covers[i] = first ? c : 4;
                first = false;
              }
            }

            const auto area_mod =
                arc_area(current_arcs.at(c).first, current_arcs.at(c).second,
                         face_normal, face_center);
            area[previous_vertices.at(p)] += area_mod;
            area[current_vertices.at(c)] -= area_mod;

            const auto new_edge = add_edge(previous_vertices[p], current_vertices[c], graph).first;
	    put(edge_level, new_edge, level);
	    put(arc_list, new_edge, ArcList{Arc(face_center, current_arcs.at(c).first, current_arcs.at(c).second, face_normal)});

            if (previous_edges.empty() ||
                previous_edges.at(p).first != current_edges.at(c).first)
              add_edge_between_faces(
                  to_halfedge[opposites.at(current_edges.at(c).first)],
                  from_halfedge[halfedges.at(current_edges.at(c).first)],
                  level - 1, previous_vertices.at(p), equivalent, graph);
            add_edge_between_faces(
                to_halfedge[opposites.at(current_edges.at(c).first)],
                from_halfedge[halfedges.at(current_edges.at(c).first)], level,
                current_vertices.at(c), equivalent, graph);
            if (previous_edges.empty() ||
                previous_edges.at(p).second != current_edges.at(c).second)
              add_edge_between_faces(
                  from_halfedge[opposites.at(current_edges.at(c).second)],
                  to_halfedge[halfedges.at(current_edges.at(c).second)],
                  level - 1, previous_vertices.at(p), equivalent, graph);
            add_edge_between_faces(
                from_halfedge[opposites.at(current_edges.at(c).second)],
                to_halfedge[halfedges.at(current_edges.at(c).second)], level,
                current_vertices.at(c), equivalent, graph);
          }
        }

        vector<Point> polygon;
        if (to_cover != 3) {
          polygon.push_back(previous_arcs.at(p).second);
          polygon.push_back(previous_arcs.at(p).first);
        }
        const int start = previous_edges.empty() ? current_edges.at(0).first
	  : previous_edges.at(p).first;
        for (const auto i : irange(0ul, to_cover)) {
          const auto vertex = (start + i) % 3;
          switch (covers[vertex]) {
          case 3:
            polygon.push_back(get(point, target(halfedges.at(vertex), mesh)));
            break;
          case 0:
          case 1:
          case 2:
            polygon.push_back(current_arcs.at(covers[vertex]).first);
            polygon.push_back(current_arcs.at(covers[vertex]).second);
            break;
          }
        }
        area[previous_vertices.at(p)] += polygon_area(polygon, face_normal, origin);

        for (const auto e : irange(0, 3))
          if (is_subset(previous_covers.at(p), edge_masks[e]) && (previous_edges.empty() || (previous_edges[p].first != e && previous_edges[p].second != e)))
            add_edge_between_faces(to_halfedge[opposites.at(e)],
                                   to_halfedge[halfedges.at(e)], level - 1,
                                   previous_vertices.at(p), equivalent, graph);
      }
      previous_arcs = move(current_arcs);
      previous_covers = move(current_covers);
      previous_edges = move(current_edges);
      previous_vertices = move(current_vertices);
    }

    for (const auto p : irange(0ul, previous_covers.size())) {
      vector<Point> polygon;
      const auto to_cover = previous_covers.at(p).count();
      if (to_cover != 3) {
        polygon.push_back(previous_arcs.at(p).second);
        polygon.push_back(previous_arcs.at(p).first);
      }
      const int start = previous_edges.empty() ? 0 : previous_edges.at(p).first;
      for (const auto i : irange(0ul, to_cover)) {
        const auto vertex = (start + i) % 3;
        if (previous_covers.at(p)[vertex])
          polygon.push_back(get(point, target(halfedges.at(vertex), mesh)));
      }
      area[previous_vertices.at(p)] += polygon_area(polygon, face_normal, origin);

      for (const auto e : irange(0, 3))
        if (is_subset(previous_covers.at(p), edge_masks[e]) && (previous_edges.empty() || (previous_edges[p].first != e && previous_edges[p].second != e)))
          add_edge_between_faces(to_halfedge[opposites.at(e)],
                                 to_halfedge[halfedges.at(e)], next_level - 1,
                                 previous_vertices.at(p), equivalent, graph);
    }
  }
}
}; // namespace shape

#endif // INTERSECT_FACES_HPP
