
/// @author Keith Jeffery

#include "base/MemoryArena.h"
#include "math/Math.h"

#include <cstdlib>
#include <iostream>

#define UTEST_ASSERT(x)                                                        \
    if (!(x)) {                                                                \
        std::cerr << "Test failed:" __FILE__ << ':' << __LINE__ << '\n' << #x; \
        exit(EXIT_FAILURE);                                                    \
    }

#define UTEST_EQUALS(x, y)                                                     \
    if (!(x == y)) {                                                           \
        std::cerr << "Test failed:" __FILE__ << ':' << __LINE__ << '\n' << #x; \
        exit(EXIT_FAILURE);                                                    \
    }

#define UTEST_FLOAT_EQUALS(x, y)                                               \
    if (!float_compare(x, y)) {                                                \
        std::cerr << "Test failed:" __FILE__ << ':' << __LINE__ << '\n' << #x; \
        exit(EXIT_FAILURE);                                                    \
    }

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
        auto p8 = arena.arena_new<S8>(11);
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

void run_tests()
{
    test_memory_arena();
}