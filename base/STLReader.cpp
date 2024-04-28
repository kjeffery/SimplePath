
#include "Endian.h"
#include "STLReader.h"
#include "../Shapes/Triangle.h"

#include <array>
#include <string_view>
#include <fstream>
#include <map>

namespace sp {
namespace {
struct Face
{
    std::array<unsigned, 3> vertex_indices;
    Normal3                 face_normal;
};

class VertexIndexer
{
public:
    auto operator()(const Point3& p) -> std::size_t
    {
        if (const auto it = m_vertex_to_index.find(p); it != m_vertex_to_index.end()) {
            return it->second;
        }
        m_vertex_to_index.emplace(p, m_next_index);
        const auto index = m_next_index;
        ++m_next_index;
        return index;
    }

private:
    std::size_t                   m_next_index{ 0 };
    std::map<Point3, std::size_t> m_vertex_to_index;
};

auto read_ascii_stl(std::istream& ins, const AffineTransformation& object_to_world) -> Mesh
{
    ins.seekg(0);
    assert(!"Not implemented yet");
    return { std::vector<std::size_t>{}, std::vector<Point3>{}, std::vector<Normal3>{}, object_to_world };
}

auto read_binary_stl(std::istream& ins, const AffineTransformation& object_to_world) -> Mesh
{
    ins.seekg(0);

    VertexIndexer vertex_indexer;

    std::vector<Point3>      vertices;
    std::vector<std::size_t> vertex_indices;
    std::vector<Face>        faces;

    constexpr std::size_t                binary_header_size = 80;
    std::array<char, binary_header_size> header;
    ins.read(header.data(), binary_header_size);

    std::uint32_t num_triangles;
    ins.read(reinterpret_cast<char*>(std::addressof(num_triangles)), sizeof(num_triangles));
    num_triangles = little_endian_to_native(num_triangles);

    // foreach triangle                      - 50 bytes:
    // REAL32[3] – Normal vector             - 12 bytes
    // REAL32[3] – Vertex 1                  - 12 bytes
    // REAL32[3] – Vertex 2                  - 12 bytes
    // REAL32[3] – Vertex 3                  - 12 bytes
    // UINT16    – Attribute byte count      -  2 bytes
    // end

    for (std::uint32_t i = 0; i < num_triangles; ++i) {
        static_assert(sizeof(float) == 4);
        float normal_x;
        float normal_y;
        float normal_z;

        ins.read(reinterpret_cast<char*>(&normal_x), sizeof(normal_x));
        ins.read(reinterpret_cast<char*>(&normal_y), sizeof(normal_y));
        ins.read(reinterpret_cast<char*>(&normal_z), sizeof(normal_z));

        Face face;
        face.face_normal = Normal3{ normal_x, normal_y, normal_z };

        float vertex_x;
        float vertex_y;
        float vertex_z;

        for (std::size_t j = 0; j < 3; ++j) {
            ins.read(reinterpret_cast<char*>(&vertex_x), sizeof(vertex_x));
            ins.read(reinterpret_cast<char*>(&vertex_y), sizeof(vertex_y));
            ins.read(reinterpret_cast<char*>(&vertex_z), sizeof(vertex_z));

            const Point3 vertex{ vertex_x, vertex_y, vertex_z };
            const auto   index = vertex_indexer(vertex);
            if (index >= vertices.size()) {
                vertices.push_back(vertex);
            }
            face.vertex_indices[j] = index;
            vertex_indices.push_back(index);
        }

        std::uint16_t attributes;
        ins.read(reinterpret_cast<char*>(&attributes), sizeof(attributes));

        if (is_zero(face.face_normal)) {
            const Vector3 edge0 = vertices.at(face.vertex_indices[1]) - vertices.at(face.vertex_indices[0]);
            const Vector3 edge1 = vertices.at(face.vertex_indices[2]) - vertices.at(face.vertex_indices[0]);
            face.face_normal    = Normal3{ cross(edge0, edge1) };
        }
        if (is_zero(face.face_normal)) {
            LOG_INFO("Encountered zero-area face. Skipping");
            continue;
        }
        face.face_normal = normalize(face.face_normal);
        faces.push_back(face);
    }

    // Calculate vertex normals from the face normals.
    std::vector vertex_normals(vertices.size(), Normal3{ 0.0f, 0.0f, 0.0f });
    std::for_each(std::execution::unseq, faces.cbegin(), faces.cend(), [&vertex_normals](const auto& f) {
        for (std::size_t i = 0; i < 3; ++i) {
            vertex_normals.at(f.vertex_indices[i]) += f.face_normal;
        }
    });

    std::transform(std::execution::par_unseq,
                   vertex_normals.cbegin(),
                   vertex_normals.cend(),
                   vertex_normals.begin(),
                   [](const auto& n) {
                       if (n != Normal3{ 0.0f, 0.0f, 0.0f }) {
                           return normalize(n);
                       } else {
                           LOG_WARNING("Found invalid normal");
                           return Normal3{ 0.0f, 1.0f, 0.0f };
                       }
                   });

    return Mesh{ std::move(vertex_indices), std::move(vertices), std::move(vertex_normals), object_to_world };
}
} // anonymous namespace

auto read_stl(const std::filesystem::path& file_name, const AffineTransformation& object_to_world) -> Mesh
{
    using namespace std::literals;

    LOG_DEBUG("Parsing STL ", file_name);

    if (!exists(file_name)) {
        throw std::runtime_error(std::format("File \"{}\" does not exist", file_name.string())); // TODO: STLException
    }
    if (!is_regular_file(file_name)) {
        throw std::runtime_error(std::format("File \"{}\" is not a regular file", file_name.string())); // TODO: STLException
    }

    std::ifstream ins(file_name, std::ios::binary);
    ins.exceptions(std::ios_base::badbit);

    constexpr auto                        ascii_header = "solid"sv;
    std::array<char, ascii_header.size()> header;
    ins.read(header.data(), ascii_header.size());
    if (std::strncmp(header.data(), ascii_header.data(), ascii_header.size()) == 0) {
        return read_ascii_stl(ins, object_to_world);
    }
    return read_binary_stl(ins, object_to_world);
}
} // namespace sp
