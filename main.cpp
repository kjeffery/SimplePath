///@author Keith Jeffery

#include "base/FileParser.h"
#include "base/not_null.h"
#include "base/ProgressBar.h"
#include "base/Tile.h"
#include "base/TileScheduler.h"
#include "base/Util.h"
#include "Cameras/Camera.h"
#include "Image/Image.h"
#include "Integrators/Integrator.h"
#include "math/Angles.h"
#include "math/LinearSpace3x3.h"
#include "math/Sampler.h"
#include "math/Vector3.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <ranges>
#include <string>
// #include <syncstream>
#include <thread>
#include <tuple>

namespace fs = std::filesystem;

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

sp::Sampler get_pixel_sampler(std::uint32_t x, std::uint32_t y)
{
    // TODO: if num passes == 1, get set
    return sp::Sampler::create_new_sequence((x << 16u) | y);
}

// TODO: should this be a member of Scene?
void render_thread(sp::Image&            image,
                   const unsigned        num_pixel_samples,
                   const sp::Scene&      scene,
                   const sp::Integrator& integrator,
                   sp::TileScheduler&    scheduler,
                   sp::ProgressBar&      progress_bar)
{
    while (auto scheduled_tile = scheduler.get_next_tile()) {
        const auto& tile = scheduled_tile->tile;
        // The tile iterators iterate over the entire collection of pixels for a full tile, regardless of clipping. We
        // use a filter to skip the pixels we're not interested in.
        auto in_tile = [&tile](const sp::Point2i& p) noexcept { return contains(tile, p); };
        for (auto p : std::views::all(tile) | std::views::filter(in_tile)) {
            auto sampler = get_pixel_sampler(p.x, p.y);
            for (unsigned i = 0; i < num_pixel_samples; ++i) {
                const auto    sample = sampler.get_next_2D();
                const sp::Ray ray    = scene.m_camera->generate_ray(p.x + sample.x, p.y + sample.y);
                // std::osyncstream(std::cout) << ray << '\n';
                image(p.x, p.y) += integrator.integrate(ray, scene, sampler);
            }
            image(p.x, p.y) /= num_pixel_samples;
        }
        progress_bar.update();
        progress_bar.draw();
    }
}

void render(unsigned num_threads, unsigned num_pixel_samples, const sp::Scene& scene)
{
    constexpr int num_passes = 1;

    sp::Image image(scene.image_width, scene.image_height, sp::RGB::black());

    // sp::MandelbrotIntegrator     integrator(scene.image_width, scene.image_height);
    sp::BruteForceIntegrator     integrator;
    sp::ColumnMajorTileScheduler scheduler{ scene.image_width, scene.image_height, num_passes };
    sp::ProgressBar              progress_bar(scheduler.get_num_tiles() * num_passes);
    std::vector<std::jthread>    threads;
    threads.reserve(num_threads);

    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back(&render_thread,
                             std::ref(image),
                             num_pixel_samples,
                             std::cref(scene),
                             std::cref(integrator),
                             std::ref(scheduler),
                             std::ref(progress_bar));
    }

    for (auto& t : threads) {
        assert(t.joinable());
        t.join();
    }

    sp::write(scene.output_file_name, image);
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

void enable_pretty_printing(std::ostream& outs)
{
    // Enable pretty-printing of our types.
    outs.iword(sp::k_pretty_print_key) = 1;
}

void print_usage(std::string_view exe_name)
{
    std::cout << "Usage: " << exe_name << "[--threads <n>] <filename>\n";
}

int main(const int argc, const char* const argv[])
{
    using namespace std::literals;
    enable_pretty_printing(std::cout);

    unsigned int num_threads       = std::thread::hardware_concurrency();
    unsigned int num_pixel_samples = 1u;

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
            } else if (arg == "--help"sv || arg == "-h"sv) {
                print_usage(argv[0]);
                return EXIT_SUCCESS;
            } else if (arg == "--threads"sv) {
                constexpr int num_args = 1;
                if (i + num_args >= argc) {
                    std::cerr << "Expected additional argument to '--threads'";
                    return EXIT_FAILURE;
                }
                std::tie(num_threads) = parse_args<unsigned>(argv + i + 1);
                i += num_args;
            } else if (arg == "--samples"sv) {
                constexpr int num_args = 1;
                if (i + num_args >= argc) {
                    std::cerr << "Expected additional argument to '--samples'";
                    return EXIT_FAILURE;
                }
                std::tie(num_pixel_samples) = parse_args<unsigned>(argv + i + 1);
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

    if (num_threads <= 0) {
        constexpr unsigned default_num_threads = 4u;
        std::cout << "Unable to determine thread count. Arbitrary choosing " << default_num_threads << ".\n";
        num_threads = default_num_threads;
    }

    try {
        using namespace sp::literals;
        const sp::Scene scene = parse_scene_file(file_path);
        render(num_threads, num_pixel_samples, scene);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error\n";
        return EXIT_FAILURE;
    }
}
