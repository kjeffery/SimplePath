
/// @author Keith Jeffery

#include "base/MemoryArena.h"
#include "Lights/Light.h"
#include "materials/Material.h"
#include "math/AffineSpace.h"
#include "math/Math.h"

#include <cstdlib>
#include <iostream>

#include "math/Sampler.h"
#include "math/Sampler.h"

using namespace sp;

#define UTEST_ASSERT(x)                                                        \
    if (!(x)) {                                                                \
        std::cerr << "Test failed:" __FILE__ << ':' << __LINE__ << '\n' << #x; \
        assert(false);                                                         \
        exit(EXIT_FAILURE);                                                    \
    }                                                                          \
    struct dummy

#define UTEST_EQUALS(x, y)                                                     \
    if (!(x == y)) {                                                           \
        std::cerr << "Test failed:" __FILE__ << ':' << __LINE__ << '\n' << #x; \
        assert(false);                                                         \
        exit(EXIT_FAILURE);                                                    \
    }                                                                          \
    struct dummy

#define UTEST_FLOAT_EQUALS(x, y)                                               \
    if (!float_compare(x, y)) {                                                \
        std::cerr << "Test failed:" __FILE__ << ':' << __LINE__ << '\n' << #x; \
        assert(false);                                                         \
        exit(EXIT_FAILURE);                                                    \
    }                                                                          \
    struct dummy

#define UTEST_FLOAT_EPSILON(x, y, eps)                                         \
    if (!float_compare_epsilon(x, y, eps)) {                                   \
        std::cerr << "Test failed:" __FILE__ << ':' << __LINE__ << '\n' << #x; \
        assert(false);                                                         \
        exit(EXIT_FAILURE);                                                    \
    }                                                                          \
    struct dummy

void do_test_memory_arena(std::size_t size)
{
    struct alignas(8) S8
    {
        explicit S8(int q)
        : x(q)
        {
        }

        int x;
    };

    struct alignas(16) S16
    {
        explicit S16(int q)
        : x(q)
        {
        }

        int x;
    };

    struct alignas(32) S32
    {
        explicit S32(int q)
        : x(q)
        {
        }

        int x;
    };

    struct alignas(64) S64
    {
        explicit S64(int q)
        : x(q)
        {
        }

        int x;
    };

    sp::MemoryArena arena(size);

    constexpr int passes = 64;
    for (int pass = 0; pass < passes; ++pass) {
        arena.release_all();
        auto p8  = arena.arena_new<S8>(11);
        auto p16 = arena.arena_new<S16>(13);
        auto p32 = arena.arena_new<S32>(17);
        auto p64 = arena.arena_new<S64>(19);
        UTEST_ASSERT(sp::is_aligned(p8, 8));
        UTEST_ASSERT(sp::is_aligned(p16, 16));
        UTEST_ASSERT(sp::is_aligned(p32, 32));
        UTEST_ASSERT(sp::is_aligned(p64, 64));
        UTEST_EQUALS(p8->x, 11);
        UTEST_EQUALS(p16->x, 13);
        UTEST_EQUALS(p32->x, 17);
        UTEST_EQUALS(p64->x, 19);
    }
}

void test_memory_arena()
{
    do_test_memory_arena(1UL);
    do_test_memory_arena(2UL);
    do_test_memory_arena(4UL);
    do_test_memory_arena(8UL);
    do_test_memory_arena(16UL);
    do_test_memory_arena(32UL);
    do_test_memory_arena(64UL);
    do_test_memory_arena(128UL);
    do_test_memory_arena(1024UL);
    do_test_memory_arena(2048UL);
    do_test_memory_arena(4096UL);
    do_test_memory_arena(8192UL);
}

void do_test_material(const sp::Material& material, const Normal3& normal)
{
    constexpr int n_samples = 1024;
    MemoryArena   arena;
    auto          sampler = IncoherentSampler::create_new_sequence(sp::Seed{ 999 });

    ONB onb = ONB::from_v(normal);

    int valid_samples = 0;

    float pdf_sum = 0.0f;
    for (int sn = 0; sn < n_samples; ++sn) {
        arena.release_all();
        const Vector3 wo     = onb.to_world(sample_to_uniform_hemisphere(sampler.get_next_2D()));
        const auto    result = material.sample(arena, wo, normal, sampler);

        if (result.pdf > 0.0f && result.color != RGB::black()) {
            const auto pdf = material.pdf(arena, wo, result.direction, normal, sampler);
            pdf_sum += pdf;
            const auto color = material.eval(arena, wo, result.direction, normal, sampler);
            UTEST_FLOAT_EPSILON(pdf, result.pdf, 0.1f);
            UTEST_ASSERT(compare_epsilon(color, result.color, 0.1f));
            ++valid_samples;
        }
    }

    const float valid_ratio = static_cast<float>(valid_samples) / static_cast<float>(n_samples);
    LOG_DEBUG("Valid sample percentage: ", valid_ratio * 100.0f);
    // UTEST_FLOAT_EQUALS((pdf_sum/n_samples) * uniform_hemisphere_pdf(), 1.0f);
}

void test_lambertian_bxdf()
{
    const auto normal = normalize(Normal3{ 1.0f, -1.0f, 1.0f });

    OneSampleMaterial::BxDFContainer bxdfs;
    bxdfs.emplace_back(new LambertianBRDF{ sp::RGB{ 0.7f, 0.6f, 0.5f } });
    OneSampleMaterial material{ std::move(bxdfs) };
    do_test_material(material, normal);
}

void test_beckmann_bxdf(float roughness)
{
    const auto      normal = normalize(Normal3{ 1.0f, -1.0f, 1.0f });
    constexpr float ior    = 1.5f;

    OneSampleMaterial::BxDFContainer        bxdfs;
    std::unique_ptr<MicrofacetDistribution> microfacet(new BeckmannDistribution{ roughness });
    bxdfs.emplace_back(new MicrofacetReflection{ RGB::white(), std::move(microfacet), ior });
    OneSampleMaterial material{ std::move(bxdfs) };
    do_test_material(material, normal);
}

void test_glossy_material(float roughness)
{
    const auto normal   = normalize(Normal3{ 1.0f, -1.0f, 1.0f });
    const auto material = create_beckmann_glossy_material(RGB{ 0.7f, 0.6f, 0.5f }, roughness, 1.5f);
    do_test_material(material, normal);
}

void test_sphere_light()
{
    auto tr = AffineSpace::translate(Vector3{ +0.0f, +3.0f, +0.0f }) * AffineSpace::scale(Vector3{ 0.1f, 0.1f, 0.1f });
    auto ir = AffineSpace::translate(Vector3{ -0.0f, -3.0f, -0.0f }) *
            AffineSpace::scale(Vector3{ 1.0f / 0.1f, 1.0f / 0.1f, 1.0f / 0.1f });

    // auto tr = AffineSpace::translate(Vector3{ +0.0f, +3.0f, +0.0f });
    // auto ir = AffineSpace::translate(Vector3{ -0.0f, -3.0f, -0.0f });

    // auto tr = AffineSpace::scale(Vector3{ 0.1f, 0.1f, 0.1f });
    // auto ir = AffineSpace::scale(Vector3{ 1.0f / 0.1f, 1.0f / 0.1f, 1.0f / 0.1f });

    SphereLight light(RGB{ 10.0f, 10.0f, 10.0f }, tr, ir);

    const Point3 p{ -3.0f, -1.0f, 2.0f };

    constexpr int num_samples = 128;
    auto          sampler     = IncoherentSampler::create_new_set(sp::Seed{ 0 }, num_samples);
    for (int i = 0; i < num_samples; ++i) {
        const auto result = light.sample(p, Normal3{ 0.0f, 1.0f, 0.0f }, sampler.get_next_2D());
        UTEST_ASSERT(light.intersect_p(result.m_tester.m_ray, RayLimits{}));
    }
}

void test_materials()
{
    test_lambertian_bxdf();

    test_beckmann_bxdf(0.1f);
    test_beckmann_bxdf(0.2f);
    test_beckmann_bxdf(0.3f);
    test_beckmann_bxdf(0.4f);
    test_beckmann_bxdf(0.5f);
    test_beckmann_bxdf(0.6f);
    test_beckmann_bxdf(0.7f);
    test_beckmann_bxdf(0.8f);
    test_beckmann_bxdf(0.9f);
    test_beckmann_bxdf(1.0f);

    //    test_glossy_material(0.1f);
    //    test_glossy_material(0.2f);
    //    test_glossy_material(0.3f);
    //    test_glossy_material(0.4f);
    //    test_glossy_material(0.5f);
    //    test_glossy_material(0.6f);
    //    test_glossy_material(0.7f);
    //    test_glossy_material(0.8f);
    //    test_glossy_material(0.9f);
    //    test_glossy_material(1.0f);
}

void run_tests()
{
    test_memory_arena();
    test_materials();
    test_sphere_light();
}
