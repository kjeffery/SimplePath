#pragma once

///@author Keith Jeffery

#include "Scene.h"

#include <iosfwd>
#include <stdexcept>

namespace sp {
class ParsingException : public std::runtime_error
{
    static std::string create_parse_message(const std::string& what_arg, int line_number);

public:
    ParsingException(const std::string& what_arg);
    ParsingException(const std::string& what_arg, int line_number);
};

Scene parse_file(std::istream& ins);

} // namespace sp
