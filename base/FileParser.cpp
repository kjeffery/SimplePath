///@author Keith Jeffery

#include "FileParser.h"

#include "Scene.h"
#include "Util.h"

#include "../math/Vector3.h"
#include "../shapes/Plane.h"
#include "../shapes/Primitive.h"
#include "../shapes/Sphere.h"

#include <algorithm>
#include <filesystem>
#include <iostream> // TODO: temp
#include <istream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std::literals;

namespace sp {

[[nodiscard]] std::string ParsingException::create_parse_message(const std::string& what_arg, int line_number)
{
    return what_arg + " on line " + std::to_string(line_number);
}

ParsingException::ParsingException(const std::string& what_arg)
: std::runtime_error(what_arg)
{
}

ParsingException::ParsingException(const std::string& what_arg, int line_number)
: std::runtime_error(create_parse_message(what_arg, line_number))
{
}

class InternalParsingException : public ParsingException
{
public:
    using ParsingException::ParsingException;
};

namespace {

[[nodiscard]] bool is_whitespace(char c) noexcept
{
    // This depends on the current locale, and I'm assuming something ASCII-like.
    return std::isspace(static_cast<unsigned char>(c)) != 0;
}

class Token
{
    std::string m_data;

    static bool is_valid_character(char c) noexcept
    {
        // This depends on the current locale, and I'm assuming something ASCII-like.
        return c == '_' || std::isalnum(static_cast<unsigned char>(c)) != 0;
    }

public:
    friend std::istream& operator>>(std::istream& ins, Token& token)
    {
        token.m_data.clear();
        while (is_whitespace(ins.peek())) {
            char c;
            ins.get(c);
        }
        while (is_valid_character(ins.peek())) {
            char c;
            ins.get(c);
            token.m_data.push_back(c);
        }
        return ins;
    }

    operator const std::string&() const
    {
        return m_data;
    }

    operator std::string() &&
    {
        return std::move(m_data);
    }
};

char consume_character(std::istream& ins, char expected, int line = -1)
{
    using namespace std::literals;

    char c;
    ins >> c;
    if (c != expected) {
        if (line >= 0) {
            throw ParsingException("Expected '"s + expected + "' character"s, line);
        } else {
            throw ParsingException("Expected '"s + expected + "' character"s);
        }
    }
    return c;
}

struct IntermediateSceneRepresentation
{
    struct PerspectiveCamera
    {
        // Point3  origin;
        Vector3 lookat;
        float   fov;
        float   focal_distance;
    };
};

[[nodiscard]] int parse_version(std::istream& ins, int line_number)
{
    Token token;
    ins >> token;
    if (const std::string& word = token; word != "version") {
        throw InternalParsingException("Expecting version directive");
    }

    consume_character(ins, ':', line_number);

    int version;
    ins >> version;
    return version;
}

class FileParser
{
public:
    FileParser()
    {
        m_parse_function_lookup.try_emplace("environment_light", &FileParser::parse_environment_light);
        m_parse_function_lookup.try_emplace("instance", &FileParser::parse_instance);
        m_parse_function_lookup.try_emplace("material_lambertian", &FileParser::parse_material_lambertian);
        m_parse_function_lookup.try_emplace("material_layered", &FileParser::parse_material_layered);
        m_parse_function_lookup.try_emplace("material_transmissive_dielectric",
                                            &FileParser::parse_material_transmissive_dielectric);
        m_parse_function_lookup.try_emplace("mesh", &FileParser::parse_mesh);
        m_parse_function_lookup.try_emplace("perspective_camera", &FileParser::parse_perspective_camera);
        m_parse_function_lookup.try_emplace("plane", &FileParser::parse_plane);
        m_parse_function_lookup.try_emplace("primitive", &FileParser::parse_primitive);
        m_parse_function_lookup.try_emplace("scene_parameters", &FileParser::parse_scene_parameters);
        m_parse_function_lookup.try_emplace("sphere", &FileParser::parse_sphere);
        m_parse_function_lookup.try_emplace("sphere_light", &FileParser::parse_sphere_light);

#if !defined(NDEBUG)
        for (const auto& s : valid_top_level_types) {
            assert(m_parse_function_lookup.contains(s));
        }
#endif
    }

    [[nodiscard]] Scene parse(std::istream& ins);

private:
    using LineNumberContainer = std::vector<int>;
    using StringSet           = std::set<std::string, std::less<>>;
    using ParseFunction       = void (FileParser::*)(const std::string&, const LineNumberContainer&, int);

    IntermediateSceneRepresentation parse_intermediate_scene(std::istream& ins);
    void                            parse_pass(const StringSet&, std::istream&, const LineNumberContainer&);

    void parse_environment_light(const std::string&, const LineNumberContainer&, int);
    void parse_instance(const std::string&, const LineNumberContainer&, int);
    void parse_material_lambertian(const std::string&, const LineNumberContainer&, int);
    void parse_material_layered(const std::string&, const LineNumberContainer&, int);
    void parse_material_transmissive_dielectric(const std::string&, const LineNumberContainer&, int);
    void parse_mesh(const std::string&, const LineNumberContainer&, int);
    void parse_perspective_camera(const std::string&, const LineNumberContainer&, int);
    void parse_plane(const std::string&, const LineNumberContainer&, int);
    void parse_primitive(const std::string&, const LineNumberContainer&, int);
    void parse_scene_parameters(const std::string&, const LineNumberContainer&, int);
    void parse_sphere(const std::string&, const LineNumberContainer&, int);
    void parse_sphere_light(const std::string&, const LineNumberContainer&, int);

    // TODO: This could be static
    std::map<std::string, ParseFunction, std::less<>> m_parse_function_lookup;

    // Sort these for binary search.
    // clang-format off
    //static constexpr std::string valid_top_level_types[] = {
    static constexpr std::string_view valid_top_level_types[] = {
        "environment_light"sv,
        "instance"sv,
        "material_lambertian"sv,
        "material_layered"sv,
        "material_transmissive_dielectric"sv,
        "mesh"sv,
        "perspective_camera"sv,
        "plane"sv,
        "primitive"sv,
        "scene_parameters"sv,
        "sphere"sv,
        "sphere_light"sv
    };
    // clang-format on

#if __cpp_lib_constexpr_algorithms >= 201806L
    static_assert(std::ranges::is_sorted(valid_top_level_types), "We binary search this data: it needs to be sorted");
#endif

    int                     m_image_width  = 256;
    int                     m_image_height = 256;
    std::filesystem::path   m_output_file_name;
    std::unique_ptr<Camera> m_camera;
};

Scene FileParser::parse(std::istream& ins)
{
    using namespace std::literals;

    try {
        ins.exceptions(std::ios::badbit);
        const auto intermediate = parse_intermediate_scene(ins);
    } catch (const ParsingException& e) {
        throw;
    } catch (const std::exception& e) {
        std::throw_with_nested(ParsingException("Unexpected file parsing error: "s + e.what()));
    } catch (...) {
        std::throw_with_nested(ParsingException("Unexpected file parsing error"));
    }

    // TODO: temp
#if 1
    auto sphere_shape0     = std::make_shared<Sphere>(AffineSpace::translate(Vector3{ 0.0f, 1.0f, 0.0f }),
                                                  AffineSpace::translate(Vector3{ 0.0f, -1.0f, 0.0f }));
    auto sphere_shape1     = std::make_shared<Sphere>(AffineSpace::translate(Vector3{ 2.0f, 1.0f, 0.0f }),
                                                  AffineSpace::translate(Vector3{ -2.0f, -1.0f, 0.0f }));
    auto plane_shape0      = std::make_shared<Plane>(AffineSpace::translate(Vector3{ 0.0f, 0.0f, 0.0f }),
                                                AffineSpace::translate(Vector3{ 0.0f, 0.0f, 0.0f }));
    auto sphere_material0  = std::make_shared<LambertianMaterial>(RGB{ 0.8f, 0.2f, 0.0f });
    auto sphere_material1  = std::make_shared<LambertianMaterial>(RGB{ 0.1f, 0.2f, 1.0f });
    auto plane_material0   = std::make_shared<LambertianMaterial>(RGB{ 0.6f, 0.6f, 1.0f });
    auto sphere_primitive0 = std::make_shared<GeometricPrimitive>(sphere_shape0, sphere_material0);
    auto sphere_primitive1 = std::make_shared<GeometricPrimitive>(sphere_shape1, sphere_material1);
    auto plane_primitive0  = std::make_shared<GeometricPrimitive>(plane_shape0, plane_material0);

    auto env_light = std::make_shared<EnvironmentLight>(RGB{ 1.0f, 1.0f, 1.0f });

    Scene::PrimitiveContainer geometry;
    geometry.push_back(sphere_primitive0);
    geometry.push_back(sphere_primitive1);
    geometry.push_back(plane_primitive0);
    Scene::LightContainer lights;
    lights.push_back(env_light);
    Scene scene{ geometry.cbegin(), geometry.cend(), lights.cbegin(), lights.cend() };

    if (!m_camera) {
        throw ParsingException("No camera specified");
    }
    scene.m_camera         = std::move(m_camera);
    scene.output_file_name = std::move(m_output_file_name);
    scene.image_width      = m_image_width;
    scene.image_height     = m_image_height;
    return scene;
#else
    return Scene{};
#endif
}

void FileParser::parse_pass(const StringSet& active_types, std::istream& ins, const LineNumberContainer& line_numbers)
{
    for (Token token; ins;) {
        ins >> token;
        if (ins.eof()) {
            break;
        }

        consume_character(ins, '{', line_numbers[ins.tellg()]);

        // Every top-level type in the scene description has a "{}" delineated body except for "version". We have to
        // eat this body regardless of whether we parse it in this pass or not in order to set up for the next top-level
        // type.

        // We are going to look at a subsection of the stream data. We have to remember of character offset for line
        // lookups before we read our subsection and convey this to the called functions.
        const std::string& word   = token;
        const auto         offset = ins.tellg();
        std::string        body;
        std::getline(ins, body, '}');
        if (active_types.contains(word)) {
            assert(m_parse_function_lookup.contains(word));
            auto fn = m_parse_function_lookup[word];
            (this->*fn)(body, line_numbers, offset);
        }
    }
}

void FileParser::parse_environment_light(const std::string&         body,
                                         const LineNumberContainer& line_numbers,
                                         int                        line_number_character_offset)
{
}

void FileParser::parse_instance(const std::string&         body,
                                const LineNumberContainer& line_numbers,
                                int                        line_number_character_offset)
{
}

void FileParser::parse_material_lambertian(const std::string&         body,
                                           const LineNumberContainer& line_numbers,
                                           int                        line_number_character_offset)
{
}

void FileParser::parse_material_layered(const std::string&         body,
                                        const LineNumberContainer& line_numbers,
                                        int                        line_number_character_offset)
{
}

void FileParser::parse_material_transmissive_dielectric(const std::string&         body,
                                                        const LineNumberContainer& line_numbers,
                                                        int                        line_number_character_offset)
{
}

void FileParser::parse_mesh(const std::string&         body,
                            const LineNumberContainer& line_numbers,
                            int                        line_number_character_offset)
{
}

void FileParser::parse_perspective_camera(const std::string&         body,
                                          const LineNumberContainer& line_numbers,
                                          int                        line_number_character_offset)
{
    std::istringstream ins(body);
    ins.exceptions(std::ios::badbit);

    Point3  origin{ no_init };
    Point3  look_at{ no_init };
    Vector3 up{ 0.0f, 1.0f, 0.0f };
    float   fov = 45.0f;

    for (Token token; ins;) {
        ins >> token;
        if (ins.eof()) {
            break;
        }
        consume_character(ins, ':', line_numbers[line_number_character_offset + ins.tellg()]);

        const std::string& word = token;
        if (word == "origin") {
            ins >> origin;
        } else if (word == "look_at") {
            ins >> look_at;
        } else if (word == "up") {
            ins >> up;
        } else if (word == "fov") {
            ins >> fov;
        } else {
            throw ParsingException("Unknown perspective_camera attribute: " + word,
                                   line_numbers[line_number_character_offset + ins.tellg()]);
        }
    }

    m_camera.reset(
        new PerspectiveCamera{ origin, look_at, up, Angle{ Degrees{ fov } }, m_image_width, m_image_height });
}

void FileParser::parse_plane(const std::string&         body,
                             const LineNumberContainer& line_numbers,
                             int                        line_number_character_offset)
{
}

void FileParser::parse_primitive(const std::string&         body,
                                 const LineNumberContainer& line_numbers,
                                 int                        line_number_character_offset)
{
}

void FileParser::parse_scene_parameters(const std::string&         body,
                                        const LineNumberContainer& line_numbers,
                                        int                        line_number_character_offset)
{
    std::istringstream ins(body);
    ins.exceptions(std::ios::badbit);
    for (Token token; ins;) {
        ins >> token;
        if (ins.eof()) {
            break;
        }
        consume_character(ins, ':', line_numbers[line_number_character_offset + ins.tellg()]);

        const std::string& word = token;
        if (word == "output_file_name") {
            std::string file_name;
            ins >> file_name;
            m_output_file_name = trim(file_name, '"');
        } else if (word == "width") {
            ins >> m_image_width;
        } else if (word == "height") {
            ins >> m_image_height;
        } else {
            throw ParsingException("Unknown scene_parameters attribute: " + word,
                                   line_numbers[line_number_character_offset + ins.tellg()]);
        }
    }
}

void FileParser::parse_sphere(const std::string&         body,
                              const LineNumberContainer& line_numbers,
                              int                        line_number_character_offset)

{
}

void FileParser::parse_sphere_light(const std::string&         body,
                                    const LineNumberContainer& line_numbers,
                                    int                        line_number_character_offset)
{
}

// Returns the contents of the file with blank lines and comments stripped as well as an array which tells us which line
// each character came from within the file. This structure is monotonic, and could definitely be compressed if we
// cared.
auto file_to_string(std::istream& ins)
{
    std::string      file_contents;
    std::vector<int> line_numbers; // Has an entry for each character in file_contents
    int              line_number = 0;
    for (std::string line; std::getline(ins, line);) {
        ++line_number;
        auto trimmed = trim(line);
        if (trimmed.empty() || trimmed.starts_with('#')) {
            continue;
        }

        trimmed = trim_end(trimmed, '#'); // This shouldn't create an empty string: it should have been caught before.
        assert(!trimmed.empty());
        assert(trimmed.find_first_of('#') == std::string_view::npos);
        file_contents.append(trimmed);
        file_contents.push_back(' '); // We read and discard the '\n', so we need a deliminator.
        line_numbers.insert(line_numbers.end(), trimmed.size() + 1u, line_number);
    }
    return std::make_pair(file_contents, line_numbers);
}

IntermediateSceneRepresentation FileParser::parse_intermediate_scene(std::istream& ins)
{
    // We read the entire file contents into memory because it makes our lives easier. If, for some reason, the input
    // file is too large, an easy workaround is to write the cleaned lines to a temporary file.
    const auto [file_contents, line_numbers] = file_to_string(ins);
    std::istringstream cleaned_ins(file_contents);

    int version = -1;
    try {
        version = parse_version(cleaned_ins, line_numbers[cleaned_ins.tellg()]);
    } catch (const InternalParsingException&) {
        std::throw_with_nested(ParsingException("Expects version as first directive"));
    }
    if (version != 1) {
        throw ParsingException("Unable to parse version " + std::to_string(version));
    }

    const auto post_version_offset = cleaned_ins.tellg();

    // First pass to check for invalid types.
    for (Token token; cleaned_ins;) {
        cleaned_ins >> token;
        if (cleaned_ins.eof()) {
            break;
        }

        consume_character(cleaned_ins, '{', line_numbers[cleaned_ins.tellg()]);

        if (const std::string& word = token; !std::ranges::binary_search(valid_top_level_types, word)) {
            throw ParsingException("Unknown type '" + word + "'", line_numbers[cleaned_ins.tellg()]);
        }

        std::string body_unused;
        std::getline(cleaned_ins, body_unused, '}');
    }
    // First parse scene parameters
    cleaned_ins.clear();
    cleaned_ins.seekg(post_version_offset);
    // clang-format off
    const StringSet pass_types0 = {
        "scene_parameters"
    };
    // clang-format on
    parse_pass(pass_types0, cleaned_ins, line_numbers);

    // Parse non-layered materials, lights, and cameras.
    cleaned_ins.clear();
    cleaned_ins.seekg(post_version_offset);
    // clang-format off
    const StringSet pass_types1 = {
        "environment_light",
        "material_lambertian",
        "material_transmissive_dielectric",
        "mesh",
        "perspective_camera",
        "plane",
        "sphere",
        "sphere_light"
    };
    // clang-format on
    parse_pass(pass_types1, cleaned_ins, line_numbers);

    // Parse layered materials (after "basic" materials have been parsed).
    cleaned_ins.clear();
    cleaned_ins.seekg(post_version_offset);
    // clang-format off
    const StringSet pass_types2 = {
        "material_layered"
    };
    // clang-format on
    parse_pass(pass_types2, cleaned_ins, line_numbers);

    // Parse primitives and instances (after geometry has already been parsed).
    cleaned_ins.clear();
    cleaned_ins.seekg(post_version_offset);
    // clang-format off
    const StringSet pass_types3 = {
        "instance",
        "primitive"
    };
    // clang-format on
    parse_pass(pass_types3, cleaned_ins, line_numbers);

    return IntermediateSceneRepresentation{};
}

} // namespace

Scene parse_file(std::istream& ins)
{
    FileParser parser;
    return parser.parse(ins);
}

} // namespace sp