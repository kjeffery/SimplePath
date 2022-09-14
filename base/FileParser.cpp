///@author Keith Jeffery

#include "Scene.h"
#include "Util.h"

#include "../math/Vector3.h"

#include <iostream> // TODO: temp
#include <istream>
#include <sstream>
#include <string>
#include <vector>

namespace sp {

namespace {
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
    std::string word;
    ins >> word;
    if (word != "version:") {
        throw 3;
    }

    int version;
    ins >> version;
    return version;
}

IntermediateSceneRepresentation::PerspectiveCamera parse_perspective_camera(const std::string& body)
{
    std::istringstream ins(body);
    for (std::string line; std::getline(ins, line);) {
        auto trimmed = sp::trim(line);
        if (trimmed.empty() || trimmed.starts_with('#')) {
            continue;
        }

        if (trimmed.starts_with("origin")) {
            const auto start = trimmed.find_first_of(':');
            if (start == std::string_view::npos) {
                throw 3;
            }
            std::cout << trimmed.substr(start + 1) << '\n';
        }
    }
}

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
    State state{State::read_version};

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
    }

    std::istringstream cleaned_ins(file_contents);
    const auto version = parse_version(cleaned_ins);
    if (version != 1) {
        throw 4;
    }

    for (std::string word; cleaned_ins;) {
        cleaned_ins >> word;
        const auto trimmed = sp::trim_end(word, ':');

        if (trimmed.starts_with("perspective_camera")) {
            //parse_perspective_camera(body);
            //std::cout << body << '\n';
        } else if (trimmed.starts_with("material_transmissive_dielectric")) {
        } else if (trimmed.starts_with("material_lambertian")) {
        } else if (trimmed.starts_with("material_layered")) {
        } else if (trimmed.starts_with("mesh")) {
        } else if (trimmed.starts_with("sphere")) {
        } else if (trimmed.starts_with("primitive")) {
        }
    }

    return IntermediateSceneRepresentation{};
}

} // namespace

Scene parse_file(std::istream& ins)
{
    ins.exceptions(std::ios::badbit);
    const auto intermediate = parse_intermediate_scene(ins);
    return Scene{};
}

} // namespace sp