
/// @author Keith Jeffery

#include "../shapes/Triangle.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <variant>

namespace sp {

struct Face
{
    std::array<unsigned, 3> v;
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

using ReadType =
    std::variant<std::int8_t, std::uint8_t, std::int16_t, std::uint16_t, std::int32_t, std::uint32_t, float, double>;

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

template <typename T>
std::unique_ptr<TypeReader> create_reader(FileType file_type)
{
    switch (file_type) {
    case FileType::ascii:
        return std::make_unique<AsciiTypeReader<T>>();
    }
    return std::make_unique<AsciiTypeReader<T>>();
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
        f.v[i]                  = vertex_index;
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
    std::ifstream ins(file_name);
    ins.exceptions(std::ios_base::badbit);

    // Read magic "ply"
    const std::string header = read_next(ins);
    if (header != "ply") {
        throw std::runtime_error("Invalid PLY header"); // TODO: PlyException
    }

    enum class DataType
    {
        INT8,
        UINT8,
        INT16,
        UINT16,
        INT32,
        UINT32,
        FLOAT,
        DOUBLE
    };

    FileType file_type;

    // Read format:
    //     format ascii 1.0
    //     format binary_little_endian 1.0
    //     format binary_big_endian 1.0
    const std::string format = read_next(ins);
    if (format == "format ascii 1.0") {
        file_type = FileType ::ascii;
    } else if (format == "format binary_little_endian 1.0") {
        file_type = FileType ::binary_little_endian;
    } else if (format == "format binary_big_endian 1.0") {
        file_type = FileType ::binary_big_endian;
    } else {
        throw std::runtime_error("Invalid PLY format"); // TODO: PlyException
    }

    std::uint32_t num_vertices{ 0 };
    std::uint32_t num_faces{ 0 };

    std::unique_ptr<TypeReader> vertex_type_reader;
    std::unique_ptr<TypeReader> face_count_type_reader;
    std::unique_ptr<TypeReader> face_index_type_reader;

    std::string line = read_next(ins);
    while (!line.empty() && line != "end_header") {
        if (line.starts_with("element")) {
            std::istringstream line_stream(line);

            std::string element_identifier;
            line_stream >> element_identifier;
            assert(element_identifier == "element");

            std::string element_type;
            line_stream >> element_type;
            if (element_type == "vertex") {
                line_stream >> num_vertices;
                line = read_next(ins);
                while (line.starts_with("property")) {
                    if (line.ends_with(" x") || line.ends_with(" y") || line.ends_with(" z")) {
                        if (line.contains(" float ")) {
                            vertex_type_reader = create_reader<float>(file_type);
                        } else if (line.contains(" double ")) {
                            vertex_type_reader = create_reader<double>(file_type);
                        } else {
                            throw std::runtime_error("Don't know how to handle vertex type"); // TODO: PlyException
                        }
                    }
                    line = read_next(ins);
                }
                // Need to read properties until they are no more...
                // But we need to put back when not a propery.
            } else if (element_type == "face") {
                line_stream >> num_faces;
            } else {
                std::cerr << "Unknown element '" << element_type << "'. Skipping.\n";
            }
        }

        line = read_next(ins);
    }

    // Look for "element face" and "element vertex" -- Warn on any other elements
    // Look for "vertex_index" or "vertex_indices"

    return Mesh{};
}

} // namespace sp
