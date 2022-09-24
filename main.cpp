///@author Keith Jeffery

#include "base/FileParser.h"
#include "base/not_null.h"
#include "base/Tile.h"
#include "base/TileScheduler.h"
#include "base/Util.h"
#include "Image/Image.h"
#include "Integrators/Integrator.h"
#include "math/LinearSpace3x3.h"
#include "math/Vector3.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <ranges>
#include <string>
#include <thread>
#include <tuple>

namespace fs = std::filesystem;

void enable_pretty_printing(std::ostream& outs)
{
    // Enable pretty-printing of our types.
    outs.iword(sp::k_pretty_print_key) = 1;
}

void print_usage(std::string_view exe_name)
{
    std::cout << "Usage: " << exe_name << "[--threads <n>] <filename>\n";
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

void render_thread(sp::Image&                    image,
                   const sp::Scene&              scene,
                   const sp::Integrator&         integrator,
                   sp::ColumnMajorTileScheduler& scheduler)
{
    while (auto scheduled_tile = scheduler.get_next_tile()) {
        const auto& tile = scheduled_tile->tile;
        // The tile iterators iterate over the entire collection of pixels for a full tile, regardless of clipping. We
        // use a filter to skip the pixels we're not interested in.
        auto in_tile = [&tile](const sp::Point2i& p) noexcept { return contains(tile, p); };
        for (auto p : std::views::all(tile) | std::views::filter(in_tile)) {
            const sp::Ray ray{ sp::Point3(p.x, p.y, 0.0f), sp::Vector3{ 0.0f, 0.0f, -1.0f } };
            // TODO: super-sampling
            image(p.x, p.y) = integrator.integrate(ray);
        }
    }
}

void render(unsigned num_threads, const sp::Scene& scene)
{
    constexpr int num_passes   = 1;
    constexpr int image_width  = 800;
    constexpr int image_height = 600;

    sp::Image image(image_width, image_height, sp::RGB::black());

    sp::MandelbrotIntegrator     integrator(image_width, image_height);
    sp::ColumnMajorTileScheduler scheduler{ image_width, image_height, num_passes };
    std::vector<std::jthread>    threads;
    threads.reserve(num_threads);

    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back(&render_thread,
                             std::ref(image),
                             std::cref(scene),
                             std::cref(integrator),
                             std::ref(scheduler));
    }

    for (auto& t : threads) {
        assert(t.joinable());
        t.join();
    }

    sp::write("image.pfm", image);
}

template <typename First, typename... Rest>
std::tuple<First, Rest...> parse_args(const char* const argv[])
{
    return std::tuple_cat(parse_args<First>(argv), parse_args<Rest...>(argv + 1));
}

template <>
std::tuple<int> parse_args<int>(const char* const argv[])
{
    const unsigned int a = std::stoi(argv[0]);
    return std::make_tuple(a);
}

template <>
std::tuple<unsigned> parse_args<unsigned>(const char* const argv[])
{
    const unsigned int a = std::stoul(argv[0]);
    return std::make_tuple(a);
}

int main(const int argc, const char* const argv[])
{
    using namespace std::literals;
    enable_pretty_printing(std::cout);

    unsigned int num_threads = std::thread::hardware_concurrency();

    if (argc == 1) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    std::string file_path;
    try {
        for (int i = 1; i < argc; ++i) {
            std::string_view arg(argv[i]);
            if (!arg.starts_with("--")) {
                file_path = arg;
            }
            if (arg == "--threads"sv) {
                constexpr int num_args = 1;
                if (i + num_args >= argc) {
                    std::cerr << "Expected additional argument to '--threads'";
                    return EXIT_FAILURE;
                }
                std::tie(num_threads) = parse_args<unsigned>(argv + i + 1);
                i += num_args;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        print_usage(argv[0]);
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error\n";
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    try {
        const sp::Scene scene = parse_scene_file(file_path);
        render(num_threads, scene);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error\n";
        return EXIT_FAILURE;
    }
}
