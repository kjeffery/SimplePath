#pragma once

/// @author Keith Jeffery

#include <cassert>
#include <type_traits>
#include <utility>

// Based off of Microsoft's GSL not_null

namespace sp {
template <typename T>
concept is_comparable_to_nullptr = std::is_convertible_v<decltype(std::declval<T>() != nullptr), bool>;

template <typename T>
// clang-format off
requires is_comparable_to_nullptr<T> && ! std::is_same_v<std::nullptr_t, T>
class not_null
// clang-format on
{
public:
    template <typename U>
    requires std::is_convertible_v<U, T>
    constexpr not_null(U&& u)
    : m_ptr(std::forward<U>(u))
    {
        assert(m_ptr != nullptr);
    }

    constexpr not_null(T u)
    : m_ptr(std::move(u))
    {
        assert(m_ptr != nullptr);
    }

    template <typename U>
    requires std::convertible_to<U, T>
    constexpr not_null(const not_null<U>& other)
    : not_null(other.get())
    {
        assert(m_ptr != nullptr);
    }

    not_null(const not_null&)            = default;
    not_null& operator=(const not_null&) = default;

    constexpr std::conditional_t<std::is_copy_constructible_v<T>, T, const T&> get() const
    {
        assert(m_ptr != nullptr);
        return m_ptr;
    }

    constexpr operator T() const
    {
        return get();
    }

    constexpr decltype(auto) operator->() const
    {
        return get();
    }

    constexpr decltype(auto) operator*() const
    {
        return *get();
    }

    not_null(std::nullptr_t)            = delete;
    not_null& operator=(std::nullptr_t) = delete;

    // Unwanted operators: points only to a single object.
    not_null& operator++()                     = delete;
    not_null& operator--()                     = delete;
    not_null  operator++(int)                  = delete;
    not_null  operator--(int)                  = delete;
    not_null& operator+=(std::ptrdiff_t)       = delete;
    not_null& operator-=(std::ptrdiff_t)       = delete;
    void      operator[](std::ptrdiff_t) const = delete;

private:
    T m_ptr;
};

} // namespace sp
