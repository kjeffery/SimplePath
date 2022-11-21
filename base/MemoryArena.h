#pragma once

/// @author Keith Jeffery

#include "math.h"

#include <cassert>
#include <cstddef>
#include <forward_list>
#include <type_traits>

namespace sp {

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
    MemoryArena()
    : m_allocated{}
    , m_free{}
    {
        m_allocated.emplace_front(s_min_block_size);
    }

    ~MemoryArena()
    {
    }

    template <typename T, typename... Args>
    [[nodiscard]] T* allocate(Args&&... args)
    requires(std::is_trivially_destructible_v<T>)
    {
        return raw_allocate<T>(std::forward<Args>(args)...);
    }

    template <typename T, typename... Args>
    [[nodiscard]] DestroyingPointer<T> allocate(Args&&... args)
    requires(!std::is_trivially_destructible_v<T>)
    {
        return DestroyingPointer<T>{ raw_allocate<T>(std::forward<Args>(args)...) };
    }

    template <typename T>
    [[nodiscard]] T* allocate_array(std::size_t n);

    void release_all() noexcept
    {
        m_free.splice_after(m_free.before_begin(), m_allocated);
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

    static bool is_aligned(void* p, const std::size_t alignment) noexcept
    {
        return reinterpret_cast<std::uintptr_t>(p) % alignment == 0;
    }

    // static constexpr std::size_t s_min_block_size = 4096U;
    static constexpr std::size_t s_min_block_size = 64U;

    struct MemoryBlock
    {
        explicit MemoryBlock(std::size_t size)
        : m_raw_memory{ new std::byte[size] }
        , m_size{ size }
        {
        }

        ~MemoryBlock()
        {
            delete[] m_raw_memory;
        }

        std::byte*  m_raw_memory;
        std::size_t m_size;
    };

    template <typename T, typename... Args>
    [[nodiscard]] T* raw_allocate(Args&&... args)
    {
        if (T* const p = allocate_from_block<T>(std::forward<Args>(args)...)) {
            assert(is_aligned(p, alignof(T)));
            return p;
        }
        new_block<T>();
        assert(m_block_offset == 0);
        T* const p = allocate_from_block<T>(std::forward<Args>(args)...);
        assert(p);
        assert(is_aligned(p, alignof(T)));
        return p;
    }

    template <typename T, typename... Args>
    [[nodiscard]] T* allocate_from_block(Args&&... args)
    {
        assert(!m_allocated.empty());

        char* const       aligned_mem     = align_up(active_block_start() + m_block_offset, alignof(T));
        const std::size_t block_offset    = aligned_mem - active_block_start();
        const std::size_t space_available = active_block_total_size() - block_offset;
        const std::size_t space_needed    = sizeof(T);

        if (space_available < space_needed) {
            return nullptr;
        }

        m_block_offset = block_offset + space_needed;
        return new (aligned_mem) T(std::forward<Args>(args)...);
    }

    template <typename T>
    void new_block()
    {
        // This is a little pessimistic, because there's a chance that the memory is already correctly aligned and we're
        // allocating too much memory, but with a large enough block size, I don't expect this to happen often (at all).
        const std::size_t size = std::max(sizeof(T) + alignof(T), s_min_block_size);
        if (!m_free.empty() && size <= m_free.front().m_size) {
            m_allocated.splice_after(m_allocated.before_begin(), m_free, m_free.before_begin(), m_free.begin());
        } else {
            m_allocated.emplace_front(MemoryBlock(size));
        }
        m_block_offset = 0;
    }

    [[nodiscard]] char* active_block_start() noexcept
    {
        assert(!m_allocated.empty());
        return reinterpret_cast<char*>(m_allocated.front().m_raw_memory);
    }

    [[nodiscard]] const char* active_block_start() const noexcept
    {
        assert(!m_allocated.empty());
        return reinterpret_cast<char*>(m_allocated.front().m_raw_memory);
    }

    [[nodiscard]] std::size_t active_block_total_size() const noexcept
    {
        assert(!m_allocated.empty());
        return m_allocated.front().m_size;
    }

    std::forward_list<MemoryBlock> m_allocated;
    std::forward_list<MemoryBlock> m_free;

    std::size_t m_block_offset = 0U;
};

} // namespace sp
