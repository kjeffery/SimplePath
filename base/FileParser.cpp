///@author Keith Jeffery

#include "Util.h"

#include <iostream> // TODO: temp
#include <istream>
#include <string>

namespace sp {

void parse_file(std::istream& ins)
{
    ins.exceptions(std::ios::badbit);

    for (std::string line; std::getline(ins, line);) {
        const auto trimmed = sp::trim(line);
        if (trimmed.empty() || trimmed.starts_with('#')) {
            continue;
        }
        std::cout << trimmed << '\n';
    }
}

} // namespace sp