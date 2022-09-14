///@author Keith Jeffery

#include "base/FileParser.h"
#include "base/Util.h"
#include "math/Vector3.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

int main(const int argc, const char* const argv[])
{
    using namespace std::literals;

    sp::Vector3 a{33962.035f, 41563.4f, 7706.415f};
    sp::Vector3 b{-24871.969f, -30438.8f, -5643.727f};
    const auto c = cross(a, b);
                //(1556.0276, -1257.5153, -75.1656)

    const auto d = reduce_min(a);

    if (argc == 2) {
        if (std::string_view file_name{ argv[1] }; file_name == "-"sv) {
            sp::parse_file(std::cin);
        } else {
            fs::path      file_path{ file_name };
            std::ifstream ins(file_path);
            if (!ins) {
                std::cerr << "Unable to open file " << file_name << '\n';
                return EXIT_FAILURE;
            }
            sp::parse_file(ins);
        }
    }
}
