///@author Keith Jeffery

#include "Scene.h"
#include "Util.h"

#include "../math/Vector3.h"

#include <iostream> // TODO: temp
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace sp {

class ParsingException : public std::runtime_error
{
    static std::string parse_message(const std::string& what_arg, int line_number)
    {
        return what_arg + " on line " + std::to_string(line_number);
    }

public:
    ParsingException(const std::string& what_arg)
    : std::runtime_error(what_arg)
    {
    }

    ParsingException(const std::string& what_arg, int line_number)
    : std::runtime_error(parse_message(what_arg, line_number))
    {
    }
};

namespace {

class Token
{
    std::string m_data;

    static bool is_valid_character(char c) noexcept
    {
        // This depends on the current locale, and I'm assuming something ASCII-like.
        return c == '_' || std::isalnum(static_cast<unsigned char>(c)) != 0;
    }

    static bool is_space(char c) noexcept
    {
        // This depends on the current locale, and I'm assuming something ASCII-like.
        return std::isspace(static_cast<unsigned char>(c)) != 0;
    }

public:
    friend std::istream& operator>>(std::istream& ins, Token& token)
    {
        while (is_space(ins.peek())) {
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

void parse_type(std::string_view type_id, std::istream& ins)
{
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

int parse_version(std::istream& ins)
{
    Token token;
    ins >> token;
    const std::string& word = token;
    if (word != "version") {
        throw ParsingException("Expecting version directive");
    }

    consume_character(ins, ':');

    int version;
    ins >> version;
    return version;
}

#if 0
IntermediateSceneRepresentation::PerspectiveCamera parse_perspective_camera(std::istream& ins)
{
    for (std::string word; ins;) {
        ins >> word;
        const auto trimmed = sp::trim_end(word, ':');
        if (trimmed.starts_with("origin")) {
            Vector3 v{ no_init };
            ins >> v;
        }
    }
    return IntermediateSceneRepresentation::PerspectiveCamera{};
}
#else
IntermediateSceneRepresentation::PerspectiveCamera parse_perspective_camera(const std::string& body)
{
    std::istringstream ins(body);
    for (Token token; ins;) {
        ins >> token;
        consume_character(ins, ':');

        const std::string& word = token;
        if (word == "origin") {
            Vector3 v{ no_init };
            ins >> v;
        } else if (word == "lookat") {
            Vector3 v{ no_init };
            ins >> v;
        } else if (word == "fov") {
            float fov;
            ins >> fov;
        } else if (word == "focal_distance") {
            float fd;
            ins >> fd;
        }
    }
    return IntermediateSceneRepresentation::PerspectiveCamera{};
}
#endif

enum class State
{
    default_state,
    read_version,
    read_perspective_camera,
    read_material_transmissive_dielectric,
    read_material_lambertian,
    read_material_layered,
    read_mesh,
    read_sphere,
    read_primitive
};

IntermediateSceneRepresentation parse_intermediate_scene(std::istream& ins)
{
    State state{ State::read_version };

    // We read the entire file contents into memory because it makes our lives easier. If, for some reason, the input
    // file is too large, an easy workaround is to write the cleaned lines to a temporary file.
    std::string file_contents;
    for (std::string line; std::getline(ins, line);) {
        auto trimmed = sp::trim(line);
        if (trimmed.empty() || trimmed.starts_with('#')) {
            continue;
        }

        trimmed = trim_end(trimmed, '#'); // This shouldn't create an empty string: it should have been caught before.
        assert(!trimmed.empty());
        assert(trimmed.find_first_of('#') == std::string_view::npos);
        file_contents.append(trimmed);
        file_contents.push_back(' '); // We read and discard the '\n', so we need a deliminator.
    }

    std::istringstream cleaned_ins(file_contents);
    try {
        const auto version = parse_version(cleaned_ins);
        if (version != 1) {
            throw ParsingException("Unable to parse version " + std::to_string(version));
        }
    } catch (const ParsingException&) {
        throw ParsingException("Expects version as first directive");
    }

    for (Token token; cleaned_ins;) {
        cleaned_ins >> token;

        consume_character(cleaned_ins, '{');

        const std::string& word = token;
        if (std::string body; word.starts_with("perspective_camera")) {
            std::getline(cleaned_ins, body, '}');
            parse_perspective_camera(body);
        } else if (word.starts_with("material_transmissive_dielectric")) {
        } else if (word.starts_with("material_lambertian")) {
        } else if (word.starts_with("material_layered")) {
        } else if (word.starts_with("mesh")) {
        } else if (word.starts_with("sphere")) {
        } else if (word.starts_with("primitive")) {
        }
    }

    return IntermediateSceneRepresentation{};
}

} // namespace

Scene parse_file(std::istream& ins)
{
    using namespace std::literals;

    try {
        ins.exceptions(std::ios::badbit);
        const auto intermediate = parse_intermediate_scene(ins);
        return Scene{};
    } catch (const ParsingException& e) {
        throw;
    } catch (const std::exception& e) {
        throw ParsingException("Unexpected file parsing error: "s + e.what());
    } catch (...) {
        throw ParsingException("Unexpected file parsing error");
    }
}

} // namespace sp