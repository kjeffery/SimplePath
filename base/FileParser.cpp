///@author Keith Jeffery

#include "FileParser.h"

#include "Scene.h"
#include "Util.h"

#include "../math/Vector3.h"

#include <iostream> // TODO: temp
#include <istream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace sp {

std::string ParsingException::parse_message(const std::string& what_arg, int line_number)
{
    return what_arg + " on line " + std::to_string(line_number);
}

ParsingException::ParsingException(const std::string& what_arg)
: std::runtime_error(what_arg)
{
}

ParsingException::ParsingException(const std::string& what_arg, int line_number)
: std::runtime_error(parse_message(what_arg, line_number))
{
}

class InternalParsingException : public ParsingException
{
public:
    using ParsingException::ParsingException;
};

namespace {

bool is_whitespace(char c) noexcept
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

int parse_version(std::istream& ins, int line_number)
{
    Token token;
    ins >> token;
    const std::string& word = token;
    if (word != "version") {
        throw InternalParsingException("Expecting version directive");
    }

    consume_character(ins, ':', line_number);

    int version;
    ins >> version;
    return version;
}

IntermediateSceneRepresentation::PerspectiveCamera parse_perspective_camera(const std::string&      body,
                                                                            const std::vector<int>& line_numbers,
                                                                            int line_number_character_offset)
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
        if (word == "origin") {
            Vector3 v{ no_init };
            ins >> v;
        } else if (word == "look_at") {
            Vector3 v{ no_init };
            ins >> v;
        } else if (word == "fov") {
            float fov;
            ins >> fov;
        } else if (word == "focal_distance") {
            float fd;
            ins >> fd;
        } else {
            throw ParsingException("Unknown perspective_camera attribute: " + word,
                                   line_numbers[line_number_character_offset + ins.tellg()]);
        }
    }
    return IntermediateSceneRepresentation::PerspectiveCamera{};
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

IntermediateSceneRepresentation parse_intermediate_scene(std::istream& ins)
{
    // We read the entire file contents into memory because it makes our lives easier. If, for some reason, the input
    // file is too large, an easy workaround is to write the cleaned lines to a temporary file.
    const auto [file_contents, line_numbers] = file_to_string(ins);
    std::istringstream cleaned_ins(file_contents);

    int version = -1;
    try {
        version = parse_version(cleaned_ins, line_numbers[cleaned_ins.tellg()]);
    } catch (const InternalParsingException&) {
        throw ParsingException("Expects version as first directive");
    }
    if (version != 1) {
        throw ParsingException("Unable to parse version " + std::to_string(version));
    }

    for (Token token; cleaned_ins;) {
        cleaned_ins >> token;
        if (cleaned_ins.eof()) {
            break;
        }

        consume_character(cleaned_ins, '{', line_numbers[cleaned_ins.tellg()]);

        // We are going to look at a subsection of the stream data. We have to remember of character offset for line
        // lookups before we read our subsection and convey this to the called functions.
        const std::string& word = token;
        if (std::string body; word.starts_with("perspective_camera")) {
            const auto offset = cleaned_ins.tellg();
            std::getline(cleaned_ins, body, '}');
            parse_perspective_camera(body, line_numbers, offset);
        } else if (word.starts_with("material_transmissive_dielectric")) {
            std::getline(cleaned_ins, body, '}');
        } else if (word.starts_with("material_lambertian")) {
            std::getline(cleaned_ins, body, '}');
        } else if (word.starts_with("material_layered")) {
            std::getline(cleaned_ins, body, '}');
        } else if (word.starts_with("mesh")) {
            std::getline(cleaned_ins, body, '}');
        } else if (word.starts_with("sphere")) {
            std::getline(cleaned_ins, body, '}');
        } else if (word.starts_with("primitive")) {
            std::getline(cleaned_ins, body, '}');
        } else {
            throw ParsingException("Unknown type '" + word + "'", line_numbers[cleaned_ins.tellg()]);
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