
/// @author Keith Jeffery

#include "Endian.h"

#include "../shapes/Triangle.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <variant>

namespace sp {

struct Face
{
    std::array<unsigned, 3> vertex_indices;
    Normal3                 face_normal;
};

// Ply types:
// name        type        number of bytes
//---------------------------------------
// char       character                 1
// uchar      unsigned character        1
// short      short integer             2
// ushort     unsigned short integer    2
// int        integer                   4
// uint       unsigned integer          4
// float      single-precision float    4
// double     double-precision float    8

// ply
// format ascii 1.0
// comment author: Greg Turk
// comment object: another cube
// element vertex 8
// property float x
// property float y
// property float z
// property uchar red                   { start of vertex color }
// property uchar green
// property uchar blue
// element face 7
// property list uchar int vertex_index  { number of vertices for each face }
// element edge 5                        { five edges in object }
// property int vertex1                  { index to first vertex of edge }
// property int vertex2                  { index to second vertex }
// property uchar red                    { start of edge color }
// property uchar green
// property uchar blue
// end_header
// 0 0 0 255 0 0                         { start of vertex list }
// 0 0 1 255 0 0
// 0 1 1 255 0 0
// 0 1 0 255 0 0
// 1 0 0 0 0 255
// 1 0 1 0 0 255
// 1 1 1 0 0 255
// 1 1 0 0 0 255
// 3 0 1 2                           { start of face list, begin with a triangle }
// 3 0 2 3                           { another triangle }
// 4 7 6 5 4                         { now some quadrilaterals }
// 4 0 4 5 1
// 4 1 5 6 2
// 4 2 6 7 3
// 4 3 7 4 0
// 0 1 255 255 255                   { start of edge list, begin with white edge }
// 1 2 255 255 255
// 2 3 255 255 255
// 3 0 255 255 255
// 2 0 0 0 0                         { end with a single black line }

enum class FileType
{
    ascii,
    binary_little_endian,
    binary_big_endian
};

enum class DataType
{
    NOT_SET,
    INT8,
    UINT8,
    INT16,
    UINT16,
    INT32,
    UINT32,
    FLOAT,
    DOUBLE
};

using ReadType =
    std::variant<std::int8_t, std::uint8_t, std::int16_t, std::uint16_t, std::int32_t, std::uint32_t, float, double>;

DataType text_to_data_type(std::string_view s)
{
    using namespace std::literals;

    if (s == "char"sv) {
        return DataType::INT8;
    } else if (s == "int8"sv) {
        return DataType::INT8;
    } else if (s == "uchar"sv) {
        return DataType::UINT8;
    } else if (s == "uint8"sv) {
        return DataType::UINT8;
    } else if (s == "short"sv) {
        return DataType::INT16;
    } else if (s == "int16"sv) {
        return DataType::INT16;
    } else if (s == "ushort"sv) {
        return DataType::UINT16;
    } else if (s == "uint16"sv) {
        return DataType::UINT16;
    } else if (s == "int"sv) {
        return DataType::INT32;
    } else if (s == "int32"sv) {
        return DataType::INT32;
    } else if (s == "uint"sv) {
        return DataType::UINT32;
    } else if (s == "uint32"sv) {
        return DataType::UINT32;
    } else if (s == "float"sv) {
        return DataType::FLOAT;
    } else if (s == "float32"sv) {
        return DataType::FLOAT;
    } else if (s == "double"sv) {
        return DataType::DOUBLE;
    } else {
        throw std::runtime_error("Unknown data type"); // TODO: PlyException
    }
}

void check_data_consistency(DataType new_data_type, DataType old_data_type)
{
    if (old_data_type != new_data_type && old_data_type != DataType::NOT_SET) {
        throw std::runtime_error("Inconsistent data type"); // TODO: PlyException
    }
}

struct TypeReader
{
    virtual ReadType read(std::istream& ins) const = 0;
};

template <typename T>
struct AsciiTypeReader : public TypeReader
{
    ReadType read(std::istream& ins) const override
    {
        T val;
        ins >> val;
        return sp::ReadType(val);
    }
};

template <>
struct AsciiTypeReader<std::nullptr_t> : public TypeReader
{
    ReadType read(std::istream& ins) const override
    {
        return sp::ReadType();
    }
};

template <typename T>
struct LittleEndianTypeReader : public TypeReader
{
    ReadType read(std::istream& ins) const override
    {
        T val;
        ins.read(reinterpret_cast<char*>(std::addressof(val)), sizeof(val));
        return little_endian(val);
    }
};

template <>
struct LittleEndianTypeReader<std::nullptr_t> : public TypeReader
{
    ReadType read(std::istream& ins) const override
    {
        return sp::ReadType();
    }
};

template <typename T>
struct BigEndianTypeReader : public TypeReader
{
    ReadType read(std::istream& ins) const override
    {
        T val;
        ins.read(reinterpret_cast<char*>(std::addressof(val)), sizeof(val));
        return big_endian(val);
    }
};

template <>
struct BigEndianTypeReader<std::nullptr_t> : public TypeReader
{
    ReadType read(std::istream& ins) const override
    {
        return sp::ReadType();
    }
};

template <typename T>
std::unique_ptr<TypeReader> create_reader(FileType file_type)
{
    switch (file_type) {
    case FileType::ascii:
        return std::make_unique<AsciiTypeReader<T>>();
    case FileType::binary_big_endian:
        return std::make_unique<BigEndianTypeReader<T>>();
    case FileType::binary_little_endian:
        return std::make_unique<LittleEndianTypeReader<T>>();
    }
    return std::make_unique<AsciiTypeReader<T>>();
}

std::unique_ptr<TypeReader> create_reader(FileType file_type, DataType data_type)
{
    switch (data_type) {
    case DataType::NOT_SET:
        assert(!"Should not get here");
        return std::make_unique<AsciiTypeReader<nullptr_t>>();
    case DataType::INT8:
        return create_reader<std::int8_t>(file_type);
    case DataType::UINT8:
        return create_reader<std::uint8_t>(file_type);
    case DataType::INT16:
        return create_reader<std::int16_t>(file_type);
    case DataType::UINT16:
        return create_reader<std::uint16_t>(file_type);
    case DataType::INT32:
        return create_reader<std::int32_t>(file_type);
    case DataType::UINT32:
        return create_reader<std::uint32_t>(file_type);
    case DataType::FLOAT:
        return create_reader<float>(file_type);
    case DataType::DOUBLE:
        return create_reader<double>(file_type);
    default:
        assert(!"Should not get here");
        return std::make_unique<AsciiTypeReader<nullptr_t>>();
    }
    assert(!"Should not get here");
    return std::make_unique<AsciiTypeReader<nullptr_t>>();
}

Face read_face(std::istream& ins, const TypeReader& count_reader, const TypeReader& index_reader)
{
    using count_type = std::uint32_t;

    auto       count_variant = count_reader.read(ins);
    const auto vertex_count  = std::visit([](auto arg) { return static_cast<count_type>(arg); }, count_variant);
    if (vertex_count != 3) {
        std::cerr << "Encountered a non-triangular face. Skipping\n";
    }

    Face f;
    for (count_type i = 0; i < vertex_count; ++i) {
        auto       vertex_index_variant = index_reader.read(ins);
        const auto vertex_index = std::visit([](auto arg) { return static_cast<unsigned>(arg); }, vertex_index_variant);
        f.vertex_indices[i]     = vertex_index;
    }
    return f;
}

std::string read_next(std::istream& ins)
{
    for (std::string s; std::getline(ins, s);) {
        const auto line = trim(s);
        if (line.empty() || line.starts_with("comment")) {
            continue;
        }
        return std::string{ line };
    }
    return std::string{};
}

Mesh read_ply(const std::filesystem::path& file_name)
{
    using namespace std::literals;

    using ReaderPointer   = std::unique_ptr<TypeReader>;
    using AnnotatedReader = std::pair<ReaderPointer, std::string>;

    std::ifstream ins(file_name, std::ios::binary);
    ins.exceptions(std::ios_base::badbit);

    // Read magic "ply"
    const std::string header = read_next(ins);
    if (header != "ply") {
        throw std::runtime_error("Invalid PLY header"); // TODO: PlyException
    }

    FileType file_type;

    // Read format:
    //     format ascii 1.0
    //     format binary_little_endian 1.0
    //     format binary_big_endian 1.0
    const std::string format = read_next(ins);
    if (format == "format ascii 1.0") {
        file_type = FileType::ascii;
    } else if (format == "format binary_little_endian 1.0") {
        file_type = FileType::binary_little_endian;
    } else if (format == "format binary_big_endian 1.0") {
        file_type = FileType::binary_big_endian;
    } else {
        throw std::runtime_error("Invalid PLY format"); // TODO: PlyException
    }

    std::uint32_t num_vertices{ 0 };
    std::uint32_t num_faces{ 0 };

    DataType vertex_list_count_data_type = DataType::NOT_SET;
    DataType vertex_list_index_data_type = DataType::NOT_SET;

    std::vector<AnnotatedReader> readers;

    std::string line = read_next(ins);
    while (!line.empty() && line != "end_header") {
        if (line.starts_with("element")) {
            const auto element_text = split(line);

            // We're assuming that we have "element" "type" and "number"
            if (element_text.size() < 3) {
                throw std::runtime_error("Unexpected argument count to 'element'"); // TODO: PlyException
            }

            if (element_text[1] == "vertex"sv) {
                num_vertices = std::stoul(std::string(element_text[2]));
                line         = read_next(ins);
                while (line.starts_with("property")) {
                    const auto property_text = split(line);
                    if (property_text.size() == 3) {
                        const auto data_type = text_to_data_type(property_text[1]);
                        readers.emplace_back(create_reader(file_type, data_type), property_text[2]);

                        if (property_text[2] != "x"sv && property_text[2] != "y"sv && property_text[2] != "z"sv) {
                            std::cerr << "Unknown property '" << property_text[2] << "'. Skipping.\n";
                        }
                    } else {
                        std::cerr << "Unknown property format '" << line << "'. Skipping.\n";
                    }
                    line = read_next(ins);
                }

                // No more properties to read
                continue;
            } else if (element_text[1] == "face"sv) {
                num_faces = std::stoul(std::string(element_text[2]));
                line      = read_next(ins);
                while (line.starts_with("property")) {
                    const auto property_text = split(line);
                    if (property_text.size() == 1) {
                        throw std::runtime_error("Malformed face property"); // TODO: PlyException
                    }
                    if (property_text.size() == 5 && property_text[1] == "list"sv &&
                        (property_text[4] == "vertex_indices"sv || property_text[4] == "vertex_index"sv)) {
                        vertex_list_count_data_type = text_to_data_type(property_text[2]);
                        vertex_list_index_data_type = text_to_data_type(property_text[3]);
                    } else {
                        std::cerr << "Unknown property '" << line << "'. Skipping.\n";
                    }
                    line = read_next(ins);
                }

                // No more properties to read.
                continue;
            } else {
                std::cerr << "Unknown element '" << element_text[1] << "'. Skipping.\n";
            }
        }
        line = read_next(ins);
    }

    std::unique_ptr<TypeReader> vertex_count_type_reader = create_reader(file_type, vertex_list_count_data_type);
    std::unique_ptr<TypeReader> vertex_index_type_reader = create_reader(file_type, vertex_list_index_data_type);

    std::vector<Point3> vertices;
    vertices.reserve(num_vertices);

    // We're making an assumption that the vertices are listed first. I'm not sure that PLY requires this.
    // But then again, PLY is an incredibly dumb file format. It's far too flexible and could be much more constrained,
    // making it easier to read.
    for (std::uint32_t i = 0; i < num_vertices; ++i) {
        float      x;
        float      y;
        float      z;
        const auto converter = [](auto arg) { return static_cast<float>(arg); };
        for (auto& reader : readers) {
            const std::string_view name(reader.second);
            const auto             read_value = reader.first->read(ins);
            if (name == "x"sv) {
                x = std::visit(converter, read_value);
            } else if (name == "y"sv) {
                y = std::visit(converter, read_value);
            } else if (name == "z"sv) {
                z = std::visit(converter, read_value);
            }
        }
        vertices.emplace_back(x, y, z);
    }

    std::vector<std::size_t> vertex_indices;
    std::vector<Face>        faces;
    // TODO: If we end up splitting quads in the future, we may want to make this reserve twice as big.
    faces.reserve(num_faces);

    for (std::uint32_t i = 0; i < num_faces; ++i) {
        auto       vertex_count_variant = vertex_count_type_reader->read(ins);
        const auto vertex_count         = std::visit(
            [](auto arg) {
                if (!std::is_integral_v<decltype(arg)>) {
                    throw std::runtime_error("Face vertex count should be integral."); // TODO: PlyException
                }
                return static_cast<std::uint64_t>(arg);
            },
            vertex_count_variant);

        // TODO: it should be easy to support quads and split them here.
        if (vertex_count != 3) {
            std::cerr << "Encountered a non-triangular face. Skipping\n";
            for (std::uint64_t i = 0; i < vertex_count; ++i) {
                vertex_index_type_reader->read(ins);
            }
            continue;
        }
        Face f;
        for (std::uint64_t v = 0; v < vertex_count; ++v) {
            auto       vertex_index_variant = vertex_index_type_reader->read(ins);
            const auto vertex_index =
                std::visit([](auto arg) { return static_cast<unsigned>(arg); }, vertex_index_variant);
            f.vertex_indices[v] = vertex_index;
        }

        // Assumption: we've read the vertices
        // Assumption: counter-clockwise vertices
        const Vector3 edge0 = vertices.at(f.vertex_indices[1]) - vertices.at(f.vertex_indices[0]);
        const Vector3 edge1 = vertices.at(f.vertex_indices[2]) - vertices.at(f.vertex_indices[0]);
        f.face_normal       = Normal3{ cross(edge0, edge1) };
        if (sqr_length(f.face_normal) < 0.00001f) {
            std::cerr << "Encountered zero-area face. Skipping\n";
            continue;
        }
        f.face_normal = normalize(f.face_normal);
        for (std::size_t v = 0; v < 3; ++v) {
            vertex_indices.push_back(f.vertex_indices[v]);
        }
        faces.push_back(f);
    }

    // Calculate vertex normals from the face normals.
    std::vector<Normal3> vertex_normals(num_vertices, Normal3{ 0.0f, 0.0f, 0.0f });
    for (const auto& f : faces) {
        for (std::size_t i = 0; i < 3; ++i) {
            vertex_normals.at(f.vertex_indices[i]) += f.face_normal;
        }
    }

    for (auto& n : vertex_normals) {
        if (n != Normal3{ 0.0f, 0.0f, 0.0f }) {
            n = normalize(n);
        }
    }

    return Mesh{ std::move(vertex_indices), std::move(vertices), std::move(vertex_normals) };
}

} // namespace sp
