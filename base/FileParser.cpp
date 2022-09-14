///@author Keith Jeffery

#include "Scene.h"
#include "Util.h"

#include "../math/Vector3.h"

#include <iostream> // TODO: temp
#include <istream>
#include <sstream>
#include <string>

namespace sp {

namespace {
void parse_type(std::string_view type_id, std::istream& ins)
{
}

struct IntermediateSceneRepresentation
{
    struct PerspectiveCamera
    {
        //Point3  origin;
        Vector3 lookat;
        float   fov;
        float   focal_distance;
    };
};

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

IntermediateSceneRepresentation parse_intermediate_scene(std::istream& ins)
{
    for (std::string line; std::getline(ins, line);) {
        const auto trimmed = sp::trim(line);
        if (trimmed.empty() || trimmed.starts_with('#')) {
            continue;
        }

        if (std::string body; trimmed.starts_with("perspective_camera")) {
            std::getline(ins, body, '}');
            parse_perspective_camera(body);
            //std::cout << body << '\n';
        } else if (trimmed.starts_with("material_transmissive_dielectric")) {
        } else if (trimmed.starts_with("material_lambertian")) {
        } else if (trimmed.starts_with("material_layered")) {
        } else if (trimmed.starts_with("mesh")) {
        } else if (trimmed.starts_with("sphere")) {
        } else if (trimmed.starts_with("primitive")) {
        }

        //std::cout << trimmed << '\n';
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