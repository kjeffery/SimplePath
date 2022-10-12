

#pragma once

/// @author Keith Jeffery

#include "Aggregate.h"
#include "ListAccelerator.h"

#include "../math/BBox.h"

#include <iterator>
#include <vector>

namespace sp {
class BVHAccelerator : public Aggregate
{
    struct NodeBase
    {
        virtual ~NodeBase() = default;

        explicit NodeBase(BBox3 bounds) noexcept
        : m_bounds(bounds)
        {
        }

        [[nodiscard]] virtual bool
        intersect(const Ray& ray, RayLimits& limits, LightIntersection& intersection) noexcept = 0;
        [[nodiscard]] virtual bool
        intersect(const Ray& ray, RayLimits& limits, Intersection& intersection) noexcept              = 0;
        [[nodiscard]] virtual bool intersect_p(const Ray& ray, const RayLimits& limits) const noexcept = 0;

        BBox3 m_bounds;
    };

    struct NodeInternal : NodeBase
    {
        NodeInternal(BBox3 bounds, std::unique_ptr<NodeBase> left, std::unique_ptr<NodeBase> right) noexcept
        : NodeBase(bounds)
        , m_children{ std::move(left), std::move(right) }
        {
        }

        bool intersect(const Ray& ray, RayLimits& limits, LightIntersection& intersection) noexcept override
        {
            bool hit = false;

            for (int i = 0; i < 2; ++i) {
                assert(m_children[i]);
                RayLimits internal_limits{ limits };
                if (sp::intersect_p(m_children[i]->m_bounds, ray, internal_limits)) {
                    if (m_children[i]->intersect(ray, internal_limits, intersection)) {
                        hit = true;
                        // The intersection with the bounding box will modify the limit's minimum value as well, which
                        // we don't want to modify for our "real" ray limits. We should only get here if we intersect
                        // actual geometry: if we hit a bounding box and no geometry, we don't want to modify the
                        // limits.
                        limits.m_t_max = internal_limits.m_t_max;
                    }
                }
            }
            return hit;
        }

        bool intersect(const Ray& ray, RayLimits& limits, Intersection& intersection) noexcept override
        {
            bool hit = false;

            for (int i = 0; i < 2; ++i) {
                assert(m_children[i]);
                RayLimits internal_limits{ limits };
                if (sp::intersect_p(m_children[i]->m_bounds, ray, internal_limits)) {
                    if (m_children[i]->intersect(ray, internal_limits, intersection)) {
                        hit = true;
                        // The intersection with the bounding box will modify the limit's minimum value as well, which
                        // we don't want to modify for our "real" ray limits. We should only get here if we intersect
                        // actual geometry: if we hit a bounding box and no geometry, we don't want to modify the
                        // limits.
                        limits.m_t_max = internal_limits.m_t_max;
                    }
                }
            }
            return hit;
        }

        bool intersect_p(const Ray& ray, const RayLimits& limits) const noexcept override
        {
            for (int i = 0; i < 2; ++i) {
                assert(m_children[i]);
                RayLimits internal_limits{ limits };
                if (sp::intersect_p(m_children[i]->m_bounds, ray, internal_limits)) {
                    if (m_children[i]->intersect_p(ray, limits)) {
                        return true;
                    }
                }
            }
            return false;
        }

        std::unique_ptr<NodeBase> m_children[2];
    };

    struct NodeLeaf : NodeBase
    {
        template <typename Iterator>
        requires std::input_iterator<Iterator>
        NodeLeaf(BBox3 bounds, Iterator first, Iterator last)
        : NodeBase(bounds)
        , m_primitives(first, last)
        {
        }

        bool intersect(const Ray& ray, RayLimits& limits, LightIntersection& intersection) noexcept override
        {
            return m_primitives.intersect(ray, limits, intersection);
        }

        bool intersect(const Ray& ray, RayLimits& limits, Intersection& intersection) noexcept override
        {
            return m_primitives.intersect(ray, limits, intersection);
        }

        bool intersect_p(const Ray& ray, const RayLimits& limits) const noexcept override
        {
            return m_primitives.intersect_p(ray, limits);
        }

        ListAccelerator m_primitives;
    };

public:
    template <typename Iterator>
    requires std::random_access_iterator<Iterator>
    BVHAccelerator(Iterator first, Iterator last)
    : m_root{ construct(first, last) }
    {
    }

private:
    bool intersect_impl(const Ray& ray, RayLimits& limits, LightIntersection& intersection) const noexcept override
    {
        assert(m_root);
        return m_root->intersect(ray, limits, intersection);
    }

    bool intersect_impl(const Ray& ray, RayLimits& limits, Intersection& intersection) const noexcept override
    {
        assert(m_root);
        return m_root->intersect(ray, limits, intersection);
    }

    bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        assert(m_root);
        return m_root->intersect_p(ray, limits);
    }

    BBox3 get_world_bounds_impl() const noexcept override
    {
        return m_root->m_bounds;
    }

    bool is_bounded_impl() const noexcept override
    {
        return true;
    }

    // Reorders elements and returns two additional iterators (called left and right) such that:b
    // Elements in [first, left) are completely to the left of the split location
    // Elements in [right, last) are completely to the right of the split location
    // Elements in [left, right) straddle the split location
    template <typename Iterator>
    requires std::random_access_iterator<Iterator>
    auto split(std::size_t dimension, const Point3& absolute_split_location, Iterator first, Iterator last)
    {
        // Put all of the elements that are totally to the left in place
        const auto left = std::partition(first, last, [absolute_split_location, dimension](const auto& a) {
            return a->get_world_bounds().get_upper()[dimension] < absolute_split_location[dimension];
        });

        // Put all of the elements that are totally to the left in place
        const auto right = std::partition(left, last, [absolute_split_location, dimension](const auto& a) {
            return a->get_world_bounds().get_lower()[dimension] < absolute_split_location[dimension];
        });
        return std::make_pair(left, right);
    }

    template <typename Iterator>
    requires std::random_access_iterator<Iterator>
    std::unique_ptr<NodeBase> construct(Iterator first, Iterator last)
    {
        std::unique_ptr<NodeBase> result;
        const auto                num_elements = std::distance(first, last);

        BBox3 bounds;
        for (Iterator it = first; it != last; ++it) {
            assert((*it)->is_bounded());
            bounds = merge(bounds, (*it)->get_world_bounds());
        }

        if (num_elements <= k_max_leaf_elements) {
            result.reset(new NodeLeaf{ bounds, first, last });
        } else {
            // TODO: Add surface-area heuristic.
            // Create n (e.g., 16) buckets that span the bounding box in our dimension.
            // For each object, classify which buckets it spans.
            // Caculate SAH for each bucket, and split at the best.
            const auto dimension     = max_dim(bounds.size());
            const auto [left, right] = split(dimension, center(bounds), first, last);

            // Elements in [left, right) are in both sub-trees, so we do an overlapping construction here.
            auto left_child{ construct(first, right) };
            auto right_child{ construct(left, last) };
            result.reset(new NodeInternal{ bounds, std::move(left_child), std::move(right_child) });
        }
        return result;
    }

    static constexpr std::size_t k_max_leaf_elements = 4u;

    using PrimitivePointer = std::shared_ptr<const Hitable>;
    using PrimitiveList    = std::vector<PrimitivePointer>;

    std::unique_ptr<NodeBase> m_root;
};
} // namespace sp
