#pragma once

/// @author Keith Jeffery

#include <cassert>
#include <type_traits>
#include <utility>

// Based off of Microsoft's GSL not_null

namespace sp {

template< class, class = void >
struct has_pointer_member: std::false_type { };

// specialization recognizes types that do have a nested ::type member:
template< class T >
struct has_pointer_member<T, std::void_t<typename T::pointer>> : std::true_type { };

template <typename T, typename Enable = void>
struct PointerType
{
    using pointer = T*;
};

template <typename T>
struct PointerType<T, typename std::enable_if<has_pointer_member<T>::value>::type>
{
    using pointer = typename T::pointer;
};

template <typename T>
concept is_comparable_to_nullptr = std::is_convertible_v<decltype(std::declval<T>() != nullptr), bool>;

template <typename T>
// clang-format off
requires is_comparable_to_nullptr<T> && (!std::is_same_v<std::nullptr_t, T>)
class not_null
// clang-format on
{
public:
    using pointer      = typename PointerType<T>::pointer;
    using element_type = T;

    template <typename U>
    requires std::convertible_to<U, T>
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
    not_null(not_null&&)                 = default;
    not_null& operator=(not_null&&)      = default;

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
