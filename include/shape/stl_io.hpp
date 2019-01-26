#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/graph_concepts.h>
#include <array>
#include <boost/container/static_vector.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/iterator/iterator_concepts.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/range/irange.hpp>
#include <regex>
#include <vector>
#include <iterator>

#ifndef STL_IO
#define STL_IO 1

/*
 * REAL32[3] – Normal vector
 * REAL32[3] – Vertex 1
 * REAL32[3] – Vertex 2
 * REAL32[3] – Vertex 3
 * UINT16 – Attribute byte count
 */
#define STL_TRIANGLE_SIZE (4 * 3 * 4 + 2)

#define FL R"((\S+))"

namespace shape {
static const std::regex header(R"(solid.*)");
static const std::regex facet(R"(\s*facet\s+normal\s+\S+\s+\S+\s+\S+)"
                              R"(\s*outer\s+loop)"
                              R"(\s*vertex\s+)" FL R"(\s+)" FL R"(\s+)" FL
                              R"(\s*vertex\s+)" FL R"(\s+)" FL R"(\s+)" FL
                              R"(\s*vertex\s+)" FL R"(\s+)" FL R"(\s+)" FL
                              R"(\s*endloop)"
                              R"(\s*endfacet)");
static const std::regex footer(R"(\s*endsolid.*\s*)");

static float little_endian(const char *chars) {
  const auto *bytes = reinterpret_cast<const std::uint8_t *>(chars);
  const uint32_t src =
      bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
  float dst;
  memcpy(&dst, &src, sizeof(float));
  return dst;
}

template <typename Iterator, typename Mesh, typename PointMap>
bool read_STL(Iterator begin, Iterator end, Mesh &mesh, PointMap &points) {
  BOOST_CONCEPT_ASSERT(
      (boost_concepts::RandomAccessTraversalConcept<Iterator>));
  BOOST_CONCEPT_ASSERT((boost_concepts::ReadableIteratorConcept<Iterator>));
  using results = std::match_results<Iterator>;
  using iterator = std::regex_iterator<Iterator>;
  BOOST_CONCEPT_ASSERT((CGAL::concepts::MutableFaceGraph<Mesh>));
  using vertex = typename boost::graph_traits<Mesh>::vertex_descriptor;
  BOOST_CONCEPT_ASSERT((boost::WritablePropertyMapConcept<PointMap, vertex>));
  using Point = typename boost::property_traits<PointMap>::value_type;
  using Diff = typename std::iterator_traits<Iterator>::difference_type;

  using boost::irange;
  using boost::container::static_vector;
  using CGAL::Euler::add_face;
  using std::array;
  using std::distance;
  using std::forward;
  using std::invalid_argument;
  using std::map;
  using std::vector;
  using std::regex_constants::match_continuous;

  auto find_vertex = [&mesh, &points,
                      vertices = map<Point, vertex>()](auto... args) mutable {
    const Point point(args...);
    const auto vertex_iterator = vertices.lower_bound(point);
    if (vertex_iterator != vertices.end() && !(vertex_iterator->first > point))
      return vertex_iterator->second;

    const auto new_vertex = add_vertex(mesh);
    put(points, new_vertex, point);
    vertices.emplace_hint(vertex_iterator, point, new_vertex);
    return new_vertex;
  };

  results name;

  if (regex_search(begin, end, name, header, match_continuous)) {
    auto last = name[0].second;
    static_vector<vertex, 3> face_vertices;
    for (iterator it = iterator(last, end, facet, match_continuous);
         it != iterator(); ++it) {
      face_vertices.clear();
      try {
        for (const auto v : irange(1, 10, 3)) {
          face_vertices.push_back(find_vertex(
              stof((*it)[v]), stof((*it)[v + 1]), stof((*it)[v + 2])));
        }
      } catch (invalid_argument &invalid) {
        return false;
      }
      add_face(face_vertices, mesh);
      last = (*it)[0].second;
    }
    if (!regex_match(last, end, footer, match_continuous)) {
      return false;
    }
  } else if (const Diff data_size = distance(begin, end);
             data_size >= 84 && ((data_size - 84) % STL_TRIANGLE_SIZE) == 0) {
    static_vector<vertex, 3> face_vertices;
    for (const auto triangle : irange<Diff>(84l, data_size, STL_TRIANGLE_SIZE)) {
      face_vertices.clear();
      for (const auto v : irange(3, 12, 3)) {
        face_vertices.push_back(
            find_vertex(little_endian(&begin[triangle + v * 4]),
                        little_endian(&begin[triangle + (v + 1) * 4]),
                        little_endian(&begin[triangle + (v + 2) * 4])));
      }
      add_face(face_vertices, mesh);
    }
  } else {
    return false;
  }
  return true;
}
}; // namespace shape

#endif // STL_IO
