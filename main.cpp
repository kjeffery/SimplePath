///@author Keith Jeffery

#include "base/FileParser.h"
#include "base/not_null.h"
#include "base/Tile.h"
#include "base/TileScheduler.h"
#include "base/Util.h"
#include "math/LinearSpace3x3.h"
#include "math/Vector3.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <ranges>
#include <string>

namespace fs = std::filesystem;

void enable_pretty_printing(std::ostream& outs)
{
    // Enable pretty-printing of our types.
    outs.iword(sp::k_pretty_print_key) = 1;
}

void print_usage(std::string_view exe_name)
{
    std::cout << "Usage: " << exe_name << " <filename>\n";
}

sp::Scene parse_scene_file(std::string_view file_name)
{
    using namespace std::literals;

    if (file_name == "-"sv) {
        return sp::parse_file(std::cin);
    } else {
        fs::path      file_path{ file_name };
        std::ifstream ins(file_path);
        if (!ins) {
            throw sp::ParsingException{ "Unable to open file "s + file_path.string() };
        }
        return sp::parse_file(ins);
    }
}

void render()
{
    sp::ColumnMajorTileScheduler scheduler{ 9, 16, 0 };
    while (auto scheduled_tile = scheduler.get_next_tile()) {
        const auto& tile = scheduled_tile->tile;
        auto in_tile = [&tile](const sp::Point2i& p) noexcept { return contains(tile, p); };
        for (auto p : std::views::all(tile) | std::views::filter(in_tile)) {
            std::cout << p << '\n';
        }
    }
}

int main(const int argc, const char* const argv[])
{
    enable_pretty_printing(std::cout);

    if (argc == 1) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

#if 0
    sp::ColumnMajorTileScheduler scheduler{9, 16, 1};
    while (auto tile = scheduler.get_next_tile()) {
        std::cout << tile->tile.get_lower() << ' ' << tile->tile.get_upper() << '\n';
        std::cout << tile->pass << '\n';
    }

    sp::Tile tile{ sp::Point2i{ 8, 16 } };
    for (auto p : tile) {
        std::cout << p << '\n';
    }
#endif

    try {
        const sp::Scene scene = parse_scene_file(argv[1]);
        render();
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error\n";
        return EXIT_FAILURE;
    }
}
