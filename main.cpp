///@author Keith Jeffery

#include "base/FileParser.h"
#include "base/Util.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

int main(const int argc, const char* const argv[])
{
    using namespace std::literals;

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
