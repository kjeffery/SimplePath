#pragma once

/// @author Keith Jeffery

#include "../math/BBox.h"
#include "../math/Morton.h"

#include <iterator>

namespace sp {
constexpr int k_tile_dimension = 8;

static_assert(is_power_of_two(k_tile_dimension), "We assume power of two for Morton traversal");

class Tile : private BBox2i
{
public:
    using BBox2i::get_lower;
    using BBox2i::get_upper;

    explicit Tile(Point2i p) noexcept
    : BBox2i{ p, p + Point2i{ k_tile_dimension } }
    {
    }

    friend inline Tile intersect(const BBox2i& a, const Tile& b) noexcept
    {
        return Tile{ intersect(a, static_cast<BBox2i>(b)) };
    }

    friend inline Tile intersect(const Tile& a, const BBox2i& b) noexcept
    {
        return intersect(b, a);
    }

    friend bool operator==(const Tile& a, const Tile& b) noexcept = default;
    friend bool contains(const Tile& tile, const Point2i& p) noexcept
    {
        return contains_open(static_cast<BBox2i>(tile), p);
    }

private:
    explicit Tile(BBox2i bounds) noexcept
    : BBox2i(bounds)
    {
    }
};

class TilePixelIterator
{
public:
    using iterator_concept  = std::random_access_iterator_tag;
    using iterator_category = std::random_access_iterator_tag;
    using value_type        = Point2i;
    using difference_type   = int;
    using reference         = value_type&;

    TilePixelIterator() noexcept
    : m_index{ 0 }
    , m_tile{Point2i{0, 0}}
    {
    }

    explicit TilePixelIterator(Tile tile, int index = 0) noexcept
    : m_index{ index }
    , m_tile{ tile }
    {
    }

    value_type operator*() const noexcept
    {
        return dereference(m_index);
    }

    TilePixelIterator& operator++() noexcept
    {
        ++m_index;
        return *this;
    }

    TilePixelIterator operator++(int) noexcept
    {
        TilePixelIterator cp{ *this };

        this->operator++();
        return cp;
    }

    TilePixelIterator& operator--() noexcept
    {
        --m_index;
        return *this;
    }

    TilePixelIterator operator--(int) noexcept
    {
        TilePixelIterator cp{ *this };

        this->operator--();
        return cp;
    }

    TilePixelIterator& operator+=(difference_type n) noexcept
    {
        m_index += n;
        return *this;
    }

    TilePixelIterator& operator-=(difference_type n) noexcept
    {
        m_index -= n;
        return *this;
    }

    value_type operator[](difference_type n) const noexcept
    {
        return dereference(m_index + n);
    }

    friend difference_type operator-(TilePixelIterator a, TilePixelIterator b) noexcept
    {
        assert(a.m_tile == b.m_tile);
        return a.m_index - b.m_index;
    }

    friend bool operator==(const TilePixelIterator&, const TilePixelIterator&) noexcept = default;
    friend std::strong_ordering operator<=>(const TilePixelIterator& a, const TilePixelIterator& b) noexcept
    {
        assert(a.m_tile == b.m_tile);
        return a.m_index <=> b.m_index;
    }

private:
    [[nodiscard]] value_type dereference(int index) const noexcept
    {
        const auto [x, y] = morton_decode(index);
        return m_tile.get_lower() + Point2i{ x, y };
    }

    int  m_index;
    Tile m_tile;
};

inline TilePixelIterator operator+(TilePixelIterator it, int n) noexcept
{
    it += n;
    return it;
}

inline TilePixelIterator operator+(int n, TilePixelIterator it) noexcept
{
    return it + n;
}

inline TilePixelIterator operator-(TilePixelIterator it, int n) noexcept
{
    it -= n;
    return it;
}

inline TilePixelIterator begin(Tile tile) noexcept
{
    return TilePixelIterator(tile, 0);
}

inline TilePixelIterator end(Tile tile) noexcept
{
    return TilePixelIterator(tile, k_tile_dimension * k_tile_dimension);
}

inline TilePixelIterator cbegin(Tile tile) noexcept
{
    return TilePixelIterator(tile, 0);
}

inline TilePixelIterator cend(Tile tile) noexcept
{
    return TilePixelIterator(tile, k_tile_dimension * k_tile_dimension);
}
} // namespace sp