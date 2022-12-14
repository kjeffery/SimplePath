#pragma once

/// @author Keith Jeffery

#include "../math/Math.h"

#include <cassert>
#include <cstddef>
#include <forward_list>
#include <memory>
#include <type_traits>

#define USE_ARENA_ALLOCATOR 1

namespace sp {

inline bool is_aligned(void* p, const std::size_t alignment) noexcept
{
    return reinterpret_cast<std::uintptr_t>(p) % alignment == 0;
}

// Destructs, does not de-allocate.
template <typename T>
struct CallDestructor
{
    void operator()(T* t) noexcept
    {
        t->~T();
    }
};

template <typename T>
using DestroyingPointer = std::unique_ptr<T, CallDestructor<T>>;

// This class is not meant to be thread-safe. Use one per thread.
class MemoryArena
{
public:
    explicit MemoryArena(std::size_t block_size = s_default_min_block_size)
    : m_allocated{}
    , m_free{}
    , m_min_block_size(round_up_to_power_of_two(block_size))
    {
        m_allocated.emplace_front(m_min_block_size);
    }

    ~MemoryArena() = default;

    template <typename T, typename... Args>
    [[nodiscard]] T* arena_new(Args&&... args)
    requires(std::is_trivially_destructible_v<T>)
    {
        T* const mem = allocate<T>(1U);
        return new (mem) T(std::forward<Args>(args)...);
    }

    template <typename T, typename... Args>
    [[nodiscard]] DestroyingPointer<T> arena_new(Args&&... args)
    requires(!std::is_trivially_destructible_v<T>)
    {
        T* const mem = allocate<T>(1U);
        new (mem) T(std::forward<Args>(args)...);
        return DestroyingPointer<T>{ mem };
    }

    // This just allocates raw memory: there is no initialization.
    template <typename T>
    [[nodiscard]] T* allocate(std::size_t n)
    {
        if (T* const p = allocate_from_block<T>(n)) {
            assert(is_aligned(p, alignof(T)));
            return p;
        }
        new_block<T>(n);
        assert(m_block_offset == 0);
        T* const p = allocate_from_block<T>(n); // This should never fail after the call to new_block;
        assert(p);
        assert(is_aligned(p, alignof(T)));
        return p;
    }

    template <typename T>
    [[nodiscard]] T* allocate_array(std::size_t n);

    void release_all() noexcept
    {
        assert(!m_allocated.empty());
        m_free.splice_after(m_free.before_begin(), m_allocated, m_allocated.begin(), m_allocated.end());
        assert(!m_allocated.empty());
        m_block_offset = 0;
    }

private:
    static std::size_t round_up_multiple_power_two(std::size_t n, std::size_t multiple) noexcept
    {
        assert(is_power_of_two(multiple));
        const std::size_t mask   = ~(multiple - 1U);
        const std::size_t result = (n + multiple - 1U) & mask;
        assert(result % multiple == 0);
        assert(result >= n);
        assert(n + multiple > result);
        return result;
    }

    static std::byte* align_up(const std::byte* const p, const std::size_t multiple) noexcept
    {
        return reinterpret_cast<std::byte*>(round_up_multiple_power_two(reinterpret_cast<std::uintptr_t>(p), multiple));
    }

    static constexpr std::size_t s_default_min_block_size = 4096U;

    struct MemoryBlock
    {
        explicit MemoryBlock(std::size_t size)
        : m_raw_memory{ new std::byte[size] }
        , m_size{ size }
        {
        }

        MemoryBlock(MemoryBlock&& other) noexcept            = default;
        MemoryBlock(const MemoryBlock&)                      = delete;
        MemoryBlock& operator=(MemoryBlock&& other) noexcept = default;
        MemoryBlock& operator=(const MemoryBlock&)           = delete;

        std::unique_ptr<std::byte[]> m_raw_memory;
        std::size_t                 m_size;
    };

    template <typename T>
    [[nodiscard]] T* allocate_from_block(std::size_t n)
    {
        assert(!m_allocated.empty());

        std::byte* const aligned_mem = align_up(active_block_start() + m_block_offset, alignof(T));
        assert(is_aligned(aligned_mem, alignof(T)));
        assert(aligned_mem >= active_block_start());

        // Start                m_block_offset           aligned_mem             End
        // +---------------------------+--------------------+---------------------+
        // |                           |                    |                     |
        // +---------------------------+--------------------+---------------------+

        const std::size_t aligned_block_offset = aligned_mem - active_block_start();
        const std::size_t block_size           = active_block_total_size();
        if (block_size < aligned_block_offset) {
            return nullptr;
        }
        const std::size_t space_available = block_size - aligned_block_offset;
        const std::size_t space_needed    = sizeof(T) * n;

        assert(aligned_mem >= active_block_start());
        assert(aligned_block_offset >= m_block_offset);

        if (space_available < space_needed) {
            return nullptr;
        }

        m_block_offset = aligned_block_offset + space_needed;
        return reinterpret_cast<T*>(aligned_mem);
    }

    template <typename T>
    void new_block(std::size_t n)
    {
        // This is a little pessimistic, because there's a chance that the memory is already correctly aligned and we're
        // allocating too much memory, but with a large enough block size, I don't expect this to happen often (at all).
        if (const std::size_t size = std::max(sizeof(T) * n + alignof(T), m_min_block_size);
            !m_free.empty() && size <= m_free.front().m_size) {
            // Move one item from the free list to the front of the allocated list.
            // The splice_after parameters are a little unusual, in that it moves the range (first, last) (open
            // interval), so we have to bookend our node iterator.
            m_allocated.splice_after(m_allocated.before_begin(), m_free, m_free.before_begin());
        } else {
            m_allocated.emplace_front(size);
        }
        m_block_offset = 0;
    }

    [[nodiscard]] std::byte* active_block_start() noexcept
    {
        assert(!m_allocated.empty());
        return m_allocated.front().m_raw_memory.get();
    }

    [[nodiscard]] const std::byte* active_block_start() const noexcept
    {
        assert(!m_allocated.empty());
        return m_allocated.front().m_raw_memory.get();
    }

    [[nodiscard]] std::size_t active_block_total_size() const noexcept
    {
        assert(!m_allocated.empty());
        return m_allocated.front().m_size;
    }

    std::forward_list<MemoryBlock> m_allocated;
    std::forward_list<MemoryBlock> m_free;

    std::size_t m_block_offset   = 0U;
    std::size_t m_min_block_size = s_default_min_block_size;
};

#if USE_ARENA_ALLOCATOR
template <typename T>
class ArenaAllocator
{
public:
    using value_type = T;

    using is_always_equal = std::false_type;

    // I don't care if things are allocated from the same arena. They aren't deallocated anyway.
    using propagate_on_container_copy_assignment = std::false_type;
    using propagate_on_container_move_assignment = std::false_type;
    using propagate_on_container_swap            = std::false_type;

    template <typename U>
    friend class ArenaAllocator;

    explicit ArenaAllocator(MemoryArena& arena) noexcept
    : m_arena(std::addressof(arena))
    {
    }

    ArenaAllocator(ArenaAllocator&&) noexcept                 = default;
    ArenaAllocator(const ArenaAllocator&) noexcept            = default;
    ArenaAllocator& operator=(ArenaAllocator&&) noexcept      = default;
    ArenaAllocator& operator=(const ArenaAllocator&) noexcept = default;

    template <typename U>
    constexpr ArenaAllocator(const ArenaAllocator<U>& other) noexcept
    : m_arena(other.m_arena)
    {
    }

    [[nodiscard]] T* allocate(std::size_t n)
    {
        return m_arena->allocate<T>(n);
    }

    void deallocate(T*, std::size_t) const noexcept
    {
        // No-op.
    }

    friend bool operator==(const ArenaAllocator&, const ArenaAllocator&) = default;

private:
    MemoryArena* m_arena;
};

#else
template <typename T>
class ArenaAllocator
{
public:
    using value_type = T;

    using is_always_equal                        = std::true_type;
    using propagate_on_container_copy_assignment = std::false_type;
    using propagate_on_container_move_assignment = std::false_type;
    using propagate_on_container_swap            = std::false_type;

    template <typename U>
    friend class ArenaAllocator;

    explicit ArenaAllocator(MemoryArena&) noexcept
    {
    }

    template <typename U>
    constexpr ArenaAllocator(const ArenaAllocator<U>&) noexcept
    {
    }

    [[nodiscard]] T* allocate(std::size_t n)
    {
        return static_cast<T*>(::operator new(n * sizeof(T)));
    }

    void deallocate(T* p, std::size_t) noexcept
    {
        ::operator delete(p);
    }

    friend bool operator==(const ArenaAllocator&, const ArenaAllocator&) = default;

private:
    MemoryArena* m_arena;
};
#endif

} // namespace sp
