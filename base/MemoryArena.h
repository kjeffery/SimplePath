#pragma once

/// @author Keith Jeffery

#include "../math/Math.h"

#include <cassert>
#include <cstddef>
#include <forward_list>
#include <memory>
#include <type_traits>

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

    static char* align_up(const char* const p, const std::size_t multiple) noexcept
    {
        return reinterpret_cast<char*>(round_up_multiple_power_two(reinterpret_cast<std::uintptr_t>(p), multiple));
    }

    static constexpr std::size_t s_default_min_block_size = 4096U;

    struct MemoryBlock
    {
        explicit MemoryBlock(std::size_t size)
        : m_raw_memory(new std::byte[size])
        , m_size{ size }
        {
        }

        ~MemoryBlock()
        {
        }

        MemoryBlock(MemoryBlock&& other) noexcept
        : m_raw_memory(std::move(other.m_raw_memory))
        , m_size(other.m_size)
        {
            other.m_size = 0UL;
        }

        MemoryBlock(const MemoryBlock&)            = delete;
        MemoryBlock& operator=(MemoryBlock&& other) noexcept
        {
            if (this != std::addressof(other)) {
                m_raw_memory = std::move(other.m_raw_memory);
                m_size = other.m_size;
                other.m_size = 0UL;
            }
        }

        MemoryBlock& operator=(const MemoryBlock&) = delete;

        std::unique_ptr<std::byte[]> m_raw_memory;
        std::size_t                  m_size;
    };

    template <typename T>
    [[nodiscard]] T* allocate_from_block(std::size_t n)
    {
        assert(!m_allocated.empty());

        char* const aligned_mem = align_up(active_block_start() + m_block_offset, alignof(T));
        assert(is_aligned(aligned_mem, alignof(T)));

        const std::size_t block_offset    = aligned_mem - active_block_start();
        const std::size_t space_available = active_block_total_size() - block_offset;
        const std::size_t space_needed    = sizeof(T) * n;

        if (space_available < space_needed) {
            return nullptr;
        }

        m_block_offset = block_offset + space_needed;
        return reinterpret_cast<T*>(aligned_mem);
    }

    template <typename T>
    void new_block(std::size_t n)
    {
        // This is a little pessimistic, because there's a chance that the memory is already correctly aligned and we're
        // allocating too much memory, but with a large enough block size, I don't expect this to happen often (at all).
        const std::size_t size = std::max(sizeof(T) * n + alignof(T), m_min_block_size);
        if (!m_free.empty() && size <= m_free.front().m_size) {
            // Move one item from the free list to the front of the allocated list.
            // The splice_after parameters are a little unusual, in that it moves the range (first, last) (open
            // interval), so we have to bookend our node iterator.
            m_allocated.splice_after(m_allocated.before_begin(), m_free, m_free.before_begin());
        } else {
            m_allocated.emplace_front(MemoryBlock(size));
        }
        m_block_offset = 0;
    }

    [[nodiscard]] char* active_block_start() noexcept
    {
        assert(!m_allocated.empty());
        return reinterpret_cast<char*>(m_allocated.front().m_raw_memory.get());
    }

    [[nodiscard]] const char* active_block_start() const noexcept
    {
        assert(!m_allocated.empty());
        return reinterpret_cast<char*>(m_allocated.front().m_raw_memory.get());
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

template <typename T>
class ArenaAllocator
{
public:
    using value_type = T;

    using propagate_on_container_copy_assignment = std::true_type;
    using propagate_on_container_move_assignment = std::true_type;
    using propagate_on_container_swap            = std::true_type;

    template <typename U>
    friend class ArenaAllocator;

    explicit ArenaAllocator(MemoryArena& arena) noexcept
    : m_arena(std::addressof(arena))
    {
    }

    template <typename U>
    constexpr ArenaAllocator(const ArenaAllocator<U>& other) noexcept
    : m_arena(other.m_arena)
    {
    }

    [[nodiscard]] T* allocate(std::size_t n)
    {
        return m_arena->allocate<T>(n);
    }

    void deallocate(T*, std::size_t) noexcept
    {
        // No-op.
    }

    friend bool operator==(const ArenaAllocator&, const ArenaAllocator&) = default;

private:
    MemoryArena* m_arena;
};

} // namespace sp
