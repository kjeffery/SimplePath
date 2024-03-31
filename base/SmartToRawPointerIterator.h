#pragma once

/// @author Keith Jeffery

#include "not_null.h"

#include <iterator>

namespace sp::detail {
template <typename T>
decltype(auto) get_pointer(T* t) noexcept
{
    return t;
}

template <typename T>
decltype(auto) get_pointer(const T& t) noexcept
{
    return t.get();
}

template <typename T>
decltype(auto) get_pointer(const not_null<T>& t) noexcept
{
    return get_pointer(t.get());
}
} // namespace sp::detail

namespace sp {

template <typename Iterator>
requires std::input_iterator<Iterator>
class SmartToRawPointerIterator
{
    Iterator m_iterator;

public:
    using iterator_concept  = typename std::iterator_traits<Iterator>::iterator_category;
    using iterator_category = typename std::iterator_traits<Iterator>::iterator_category;
    using value_type        = typename std::iterator_traits<Iterator>::value_type::pointer;
    using difference_type   = typename std::iterator_traits<Iterator>::difference_type;

    SmartToRawPointerIterator() = default;

    explicit SmartToRawPointerIterator(Iterator it)
    : m_iterator(it)
    {
    }

    value_type operator*() const
    {
        // return m_iterator->get();
        return detail::get_pointer(*m_iterator);
    }

    SmartToRawPointerIterator& operator++()
    {
        ++m_iterator;
        return *this;
    }

    SmartToRawPointerIterator operator++(int)
    {
        SmartToRawPointerIterator cp{ *this };

        this->operator++();
        return cp;
    }

    SmartToRawPointerIterator& operator--()
    requires std::bidirectional_iterator<Iterator>
    {
        --m_iterator;
        return *this;
    }

    SmartToRawPointerIterator operator--(int)
    requires std::bidirectional_iterator<Iterator>
    {
        SmartToRawPointerIterator cp{ *this };

        this->operator--();
        return cp;
    }

    SmartToRawPointerIterator operator+=(difference_type n)
    requires std::random_access_iterator<Iterator>
    {
        m_iterator += n;
        return *this;
    }

    SmartToRawPointerIterator operator-=(difference_type n)
    requires std::random_access_iterator<Iterator>
    {
        m_iterator -= n;
        return *this;
    }

    SmartToRawPointerIterator operator[](difference_type n)
    requires std::random_access_iterator<Iterator>
    {
        return m_iterator[n];
    }

    friend difference_type operator-(SmartToRawPointerIterator a, SmartToRawPointerIterator b)
    requires std::random_access_iterator<Iterator>
    {
        return a.m_iterator - b.m_iterator;
    }

    friend bool operator==(const SmartToRawPointerIterator&, const SmartToRawPointerIterator&) = default;
};

} // namespace sp
