#pragma once

/// @author Keith Jeffery

#include "Tile.h"

#include <atomic>
#include <optional>

namespace sp {

struct ScheduledTile
{
    Tile tile;
    int  pass;
};

class TileScheduler
{
public:
    TileScheduler(int width, int height) noexcept
    : m_extents(Point2i{ 0, 0 }, Point2i{ width, height })
    {
    }

    virtual ~TileScheduler() = default;

    // This must be thread-safe.
    [[nodiscard]] std::optional<ScheduledTile> get_next_tile()
    {
        auto t = get_next_tile_impl();
        if (t) {
            t->tile = intersect(t->tile, m_extents);
        }
        return t;
    }

    template <std::size_t dim>
    int get_num_tiles() noexcept
    {
        const int s = extents<dim>(m_extents);
        return (s + (k_tile_dimension - 1)) / k_tile_dimension;
    }

    int get_num_tiles() noexcept
    {
        return get_num_tiles<0>() * get_num_tiles<1>();
    }

private:
    virtual std::optional<ScheduledTile> get_next_tile_impl() = 0;

    BBox2i m_extents;
};

class ColumnMajorTileScheduler : public TileScheduler
{
public:
    ColumnMajorTileScheduler(int width, int height, int pass_clamp) noexcept
    : TileScheduler{ width, height }
    , m_counter(0)
    , m_pass_clamp(pass_clamp)
    {
    }

private:
    std::optional<ScheduledTile> get_next_tile_impl() override
    {
        const int counter   = m_counter++;
        const int num_tiles = get_num_tiles();
        assert(num_tiles > 0);
        const int pass = counter / num_tiles;
        if (pass > m_pass_clamp) {
            return std::optional<ScheduledTile>{};
        }
        const int  index = counter % num_tiles;
        const int  w     = get_num_tiles<0>();
        const auto x     = index % w;
        const auto y     = index / w;
        return std::optional<ScheduledTile>{ std::in_place,
                                             Tile{ Point2i{ x * k_tile_dimension, y * k_tile_dimension } },
                                             pass };
    }

    std::atomic<int> m_counter;
    int              m_pass_clamp;
};

} // namespace sp
