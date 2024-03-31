///@author Keith Jeffery

#include "UnitTests.h"

#include "base/FileParser.h"
#include "base/MemoryArena.h"
#include "base/not_null.h"
#include "base/ProgressBar.h"
#include "base/Stopwatch.h"
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

namespace sp {
int k_pretty_print_key = -1;
} // namespace sp

auto create_integrator(const sp::IntegratorType type, const int image_width, const int image_height, const int min_depth,
                       const int                min_height) -> std::unique_ptr<sp::Integrator>
{
    switch (type) {
    case sp::IntegratorType::Mandelbrot: return std::make_unique<sp::MandelbrotIntegrator>(image_width, image_height);
    case sp::IntegratorType::BruteForce: return std::make_unique<sp::BruteForceIntegrator>();
    case sp::IntegratorType::BruteForceIterative: return std::make_unique<sp::BruteForceIntegratorIterative>();
    case sp::IntegratorType::BruteForceIterativeRR: return std::make_unique<sp::BruteForceIntegratorIterativeRR>();
    case sp::IntegratorType::BruteForceIterativeRRNEE: return std::make_unique<sp::BruteForceIntegratorIterativeRRNEE>();
    case sp::IntegratorType::DirectLighting: return std::make_unique<sp::DirectLightingIntegrator>();
    default: return std::make_unique<sp::BruteForceIntegratorIterative>();
    }
}

sp::Scene parse_scene_file(std::string_view file_name)
{
    using namespace std::literals;

    if (file_name == "-"sv) {
        return sp::parse_file(std::cin);
    } else {
        const fs::path      file_path{ file_name };
        std::ifstream ins(file_path);
        if (!ins) {
            throw sp::ParsingException{ "Unable to open file "s + file_path.string() };
        }
        return sp::parse_file(ins);
    }
}

auto get_pixel_sampler(std::uint32_t x, std::uint32_t y) -> sp::RSequenceSampler
{
    // TODO: if num passes == 1, get set
    return sp::RSequenceSampler::create_new_sequence(sp::Seed{ x << 16u | y });
}

auto get_integrator_sampler(std::uint32_t x, std::uint32_t y) -> sp::IncoherentSampler
{
    return sp::IncoherentSampler::create_new_sequence(sp::Seed{ (x << 16u | y) ^ 0xb0ae9d99 });
}
void render_thread(sp::Image&            image,
                   const unsigned        num_pixel_samples,
                   const sp::Scene&      scene,
                   const sp::Integrator& integrator,
                   sp::TileScheduler&    scheduler,
                   sp::ProgressBar&      progress_bar)
{
    sp::MemoryArena arena;

    while (auto scheduled_tile = scheduler.get_next_tile()) {
        const auto& tile = scheduled_tile->tile;
        // The tile iterators iterate over the entire collection of pixels for a full tile, regardless of clipping. We
        // use a filter to skip the pixels we're not interested in.
        auto in_tile = [&tile](const sp::Point2i& p) noexcept { return contains(tile, p); };
        for (auto p : std::views::all(tile) | std::views::filter(in_tile)) {
            auto pixel_sampler      = get_pixel_sampler(p.x, p.y);
            auto integrator_sampler = get_integrator_sampler(p.x, p.y);
            for (unsigned i = 0; i < num_pixel_samples; ++i) {
                arena.release_all();
                const auto       sample = pixel_sampler.get_next_2D();
                const sp::Point2 pixel_coords{ p.x + sample.x, p.y + sample.y };
                const sp::Ray    ray = scene.m_camera->generate_ray(pixel_coords.x, pixel_coords.y);
                // std::osyncstream(std::cout) << ray << '\n';
                image(p.x, p.y) += integrator.integrate(ray, scene, arena, integrator_sampler, pixel_coords);
            }
            image(p.x, p.y) /= num_pixel_samples;
        }
        progress_bar.update();
        progress_bar.draw();
    }
}

void render(const sp::Integrator& integrator, unsigned num_threads, unsigned num_pixel_samples, const sp::Scene& scene)
{
    constexpr int num_passes = 1;

    sp::Stopwatch stopwatch;

    sp::Image image(scene.image_width, scene.image_height, sp::RGB::black());

    sp::ColumnMajorTileScheduler scheduler{ scene.image_width, scene.image_height, num_passes };
    sp::ProgressBar              progress_bar(scheduler.get_num_tiles() * num_passes, "tiles");
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
    stopwatch.stop();
    std::cout << "\nElapsed time: ";
    stopwatch.print(std::cout);
}

void morton_demonstration()
{
    constexpr unsigned tile_size          = 32u;
    constexpr unsigned num_tiles_1D       = 16u;
    constexpr unsigned num_tiles          = sp::square(num_tiles_1D);
    constexpr unsigned num_pixels_1D      = tile_size * num_tiles_1D;
    constexpr unsigned frames_to_activate = 2u;
    constexpr unsigned frames_to_fade     = 50u;
    constexpr unsigned total_frames       = frames_to_activate * num_tiles + frames_to_fade;
    constexpr float    min_saturation     = 0.0f;

    constexpr sp::HSV base_color{ sp::Degrees{ 240.0f }, min_saturation, 1.0f };
    constexpr sp::HSV hit_color{ sp::Degrees{ 240.0f }, 1.0f, 1.0f };

    sp::Array2D<sp::HSV> hsv_image(num_pixels_1D, num_pixels_1D, base_color);

    auto convert_to_rgb = [](const sp::Array2D<sp::HSV>& in) {
        sp::Image img(in.width(), in.height());

        // TODO: Huh. I thought I had iterators on the Array2D class...
        // std::transform(std::execution::par_unseq, in.cbegin(), in.cend(), img.begin())
        for (sp::Image::size_type y = 0; y < in.height(); ++y) {
            for (sp::Image::size_type x = 0; x < in.width(); ++x) {
                img(x, y) = sp::to_rgb(in(x, y));
            }
        }
        return img;
    };

    auto set_tile_color = [&hsv_image, tile_size](unsigned x, unsigned y, const sp::HSV& color) {
        const auto start_x = x * tile_size;
        const auto start_y = y * tile_size;
        const auto end_x   = start_x + tile_size;
        const auto end_y   = start_y + tile_size;
        for (unsigned iy = start_y; iy < end_y; ++iy) {
            for (unsigned ix = start_x; ix < end_x; ++ix) {
                hsv_image(ix, iy) = color;
            }
        }
    };

    auto add_grid = [tile_size](sp::Image& image) {
        sp::Image result(image.width(), image.height());
        for (int i = 0; i <= image.width(); i += tile_size) {
            const int left  = std::max<int>(0, i - 1);
            const int right = std::min<int>(image.width() - 1, i);

            for (sp::Image::size_type y = 0; y < image.height(); ++y) {
                result(left, y)  = max(sp::RGB{ 0.3f }, result(left, y));
                result(right, y) = max(sp::RGB{ 0.3f }, result(right, y));
            }
        }
        for (int i = 0; i <= image.height(); i += tile_size) {
            const int top    = std::max<int>(0, i - 1);
            const int bottom = std::min<int>(image.height() - 1, i);

            for (sp::Image::size_type x = 0; x < image.width(); ++x) {
                result(x, bottom) = max(sp::RGB{ 0.3f }, result(x, bottom));
                result(x, top)    = max(sp::RGB{ 0.3f }, result(x, top));
            }
        }
        const int center_x = image.width() / 2;
        const int center_y = image.height() / 2;
        for (int x = 0; x < image.width(); ++x) {
            result(x, center_y)     = max(sp::RGB{ 0.7f }, result(x, center_y));
            result(x, center_y - 1) = max(sp::RGB{ 0.7f }, result(x, center_y - 1));
        }
        for (int y = 0; y < image.height(); ++y) {
            result(center_x, y)     = max(sp::RGB{ 0.7f }, result(center_x, y));
            result(center_x - 1, y) = max(sp::RGB{ 0.7f }, result(center_x - 1, y));
        }

        for (sp::Image::size_type y = 0; y < image.height(); ++y) {
            for (sp::Image::size_type x = 0; x < image.width(); ++x) {
                const float alpha = sp::relative_luminance(result(x, y));
                result(x, y)      = result(x, y) * alpha + image(x, y) * (1.0f - alpha);
            }
        }

        return result;
    };

    for (unsigned frame = 0; frame < total_frames; ++frame) {
        for (unsigned tile = 0; tile < num_tiles; ++tile) {
            const unsigned activation_frame = frames_to_activate * tile;
            if (frame < activation_frame) {
                continue;
            }

            const auto [x, y] = sp::morton_decode(tile);

            if (frame == activation_frame) {
                set_tile_color(x, y, hit_color);
            } else {
                assert(activation_frame < frame);

                float saturation = min_saturation;
                if (frame - activation_frame <= frames_to_fade) {
                    const unsigned frames_past_expiration = frame - activation_frame;
                    saturation                            = min_saturation + (1.0f - static_cast<float>(frames_past_expiration) /
                        static_cast<float>(frames_to_fade));
                }

                sp::HSV color = hit_color;
                color.s       = saturation;
                set_tile_color(x, y, color);
            }
        }

        // This would be nicer if GCC 10 supported format.
        sp::Image rgb_image = convert_to_rgb(hsv_image);
        rgb_image           = add_grid(rgb_image);

        std::ostringstream name_stream;
        name_stream << "morton_frames/morton_" << std::setw(4) << std::setfill('0') << frame << ".pfm";
        sp::write(name_stream.str(), rgb_image);
    }
}

template<typename First, typename... Rest>
std::tuple<First, Rest...> parse_args(const char* const argv[])
{
    return std::tuple_cat(parse_args<First>(argv), parse_args<Rest...>(argv + 1));
}

template<>
std::tuple<int> parse_args<int>(const char* const argv[])
{
    const unsigned int a = std::stoi(argv[0]);
    return std::make_tuple(a);
}

template<>
std::tuple<unsigned> parse_args<unsigned>(const char* const argv[])
{
    const unsigned int a = std::stoul(argv[0]);
    return std::make_tuple(a);
}

void pretty_print_callback(std::ios::event event, std::ios_base& b, int)
{
    if (event == std::ios::copyfmt_event) {
        b.iword(sp::k_pretty_print_key) = 1;
    }
}

void enable_pretty_printing(std::ostream& outs)
{
    // Enable pretty-printing of our types.

    // Maybe a bit of abuse of call_once? We don't care much about thread-safety here, we just want to allocate our
    // key once.
    static std::once_flag xalloc_flag;
    std::call_once(xalloc_flag, []() { sp::k_pretty_print_key = std::ios_base::xalloc(); });

    outs.iword(sp::k_pretty_print_key) = 1;
    outs.register_callback(&pretty_print_callback, sp::k_pretty_print_key);
}

void print_usage(std::string_view exe_name)
{
    std::cout << "Usage: " << exe_name << " [--threads <n>] <filename>\n";
}

int main(const int argc, const char* const argv[])
{
    using namespace std::literals;
    enable_pretty_printing(std::cout);
    enable_pretty_printing(std::cerr);

    unsigned int num_threads       = std::thread::hardware_concurrency();
    unsigned int num_pixel_samples = 1u;
    bool         run_unit_tests    = false;

    if (argc == 1) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    sp::IntegratorType integrator_type{ sp::IntegratorType::NotSpecified };

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
            } else if (arg == "--integrator"sv) {
                constexpr int num_args = 1;
                if (i + num_args >= argc) {
                    std::cerr << "Expected additional argument to '--integrator'";
                    return EXIT_FAILURE;
                }
                integrator_type = sp::string_to_integrator_type(argv[i + 1]);
                i += num_args;
            } else if (arg == "--test"sv) {
                run_unit_tests = true;
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

    if (run_unit_tests) {
        run_tests();
        return EXIT_SUCCESS;
    }

    try {
        using namespace sp::literals;
        const sp::Scene scene = parse_scene_file(file_path);

        if (integrator_type == sp::IntegratorType::NotSpecified) {
            integrator_type = scene.integrator_type;
        }
        if (integrator_type == sp::IntegratorType::NotSpecified) {
            integrator_type = sp::IntegratorType::DirectLighting;
        }

        const auto integrator = create_integrator(integrator_type, scene.image_width, scene.image_height, scene.russian_roulette_depth,
                                                  scene.max_depth);
        assert(integrator);
        render(*integrator, num_threads, num_pixel_samples, scene);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error\n";
        return EXIT_FAILURE;
    }
}
