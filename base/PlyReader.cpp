
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
    } else if (s == "uchar"sv) {
        return DataType::UINT8;
    } else if (s == "short"sv) {
        return DataType::INT16;
    } else if (s == "ushort"sv) {
        return DataType::UINT16;
    } else if (s == "int"sv) {
        return DataType::INT16;
    } else if (s == "uint"sv) {
        return DataType::UINT16;
    } else if (s == "float"sv) {
        return DataType::FLOAT;
    } else if (s == "double"sv) {
        return DataType::DOUBLE;
    } else {
        throw std::runtime_error("Unknown data type"); // TODO: PlyException
    }
}

DataType check_data_consistency(DataType new_data_type, DataType old_data_type)
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

template <typename T>
std::unique_ptr<TypeReader> create_reader(FileType file_type)
{
    switch (file_type) {
    case FileType::ascii:
        return std::make_unique<AsciiTypeReader<T>>();
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
    using namespace std::literals;

    std::ifstream ins(file_name);
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

    DataType vertex_data_type            = DataType::NOT_SET;
    DataType vertex_list_count_data_type = DataType::NOT_SET;
    DataType vertex_list_index_data_type = DataType::NOT_SET;

    std::string line = read_next(ins);
    while (!line.empty() && line != "end_header") {
        if (line.starts_with("element")) {
            const auto element_text = split(line);

            // We're assuming that we have "element" "type" and "number"
            if (element_text.size() < 3) {
                throw std::runtime_error("Unexpected argument count to 'element'"); // TODO: PlyException
            }

            if (element_text[1] == "vertex"sv) {
                num_vertices = std::strtol(element_text[2]);
                line         = read_next(ins);
                while (line.starts_with("property")) {
                    const auto property_text = split(line);
                    if (property_text.size() == 3) {
                        if (property_text[2] == "x"sv || property_text[2] == "y"sv || property_text[2] == "z"sv) {
                            const auto data_type = text_to_data_type(property_text[1]);
                            check_data_consistency(data_type, vertex_data_type);
                            vertex_data_type = data_type;
                        } else {
                            std::cerr << "Unknown property '" << line << "'. Skipping.\n";
                        }
                    } else {
                        std::cerr << "Unknown property '" << line << "'. Skipping.\n";
                    }
                }
            } else if (element_text[1] == "face"sv) {
                num_vertices = std::strtol(element_text[2]);
                line         = read_next(ins);
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
                }
            } else {
                std::cerr << "Unknown element '" << element_type << "'. Skipping.\n";
            }
        }
        line = read_next(ins);
    }
    std::unique_ptr<TypeReader> vertex_type_reader     = create_reader(file_type, vertex_data_type);
    std::unique_ptr<TypeReader> face_count_type_reader = create_reader(file_type, vertex_list_count_data_type);
    std::unique_ptr<TypeReader> face_index_type_reader = create_reader(file_type, vertex_list_index_data_type);

    // Look for "element face" and "element vertex" -- Warn on any other elements
    // Look for "vertex_index" or "vertex_indices"

    return Mesh{};
}

} // namespace sp
