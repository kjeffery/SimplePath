

#pragma once

/// @author Keith Jeffery

#include "Aggregate.h"
#include "ListAccelerator.h"

#include "../math/BBox.h"

#include <algorithm>
#include <concepts>
#include <execution>
#include <iterator>
#include <vector>

namespace sp {
class BVHAccelerator : public Aggregate
{
    struct NodeBase
    {
        virtual ~NodeBase() = default;

        explicit NodeBase(BBox3 bounds) noexcept
        : m_bounds(std::move(bounds))
        {
        }

        [[nodiscard]] virtual std::optional<LightIntersection> intersect_lights(const Ray& ray, const RayLimits& limits) noexcept = 0;
        [[nodiscard]] virtual std::optional<Intersection>      intersect(const Ray& ray, const RayLimits& limits) noexcept = 0;
        [[nodiscard]] virtual bool                             intersect_p(const Ray& ray, const RayLimits& limits) const noexcept = 0;

        BBox3 m_bounds;
    };

    struct NodeInternal : NodeBase
    {
        NodeInternal(BBox3 bounds, std::unique_ptr<NodeBase> left, std::unique_ptr<NodeBase> right) noexcept
        : NodeBase(std::move(bounds))
        , m_children{ std::move(left), std::move(right) }
        {
        }

        std::optional<LightIntersection> intersect_lights(const Ray& ray, const RayLimits& limits) noexcept override
        {
            std::optional<LightIntersection> result;

            auto intermediate_limits{ limits };
            for (const auto& child : m_children) {
                assert(child);
                if (sp::intersect_p(child->m_bounds, ray, intermediate_limits)) {
                    if (const auto intersection = child->intersect_lights(ray, intermediate_limits); intersection) {
                        intermediate_limits.m_t_max = intersection->m_distance;
                        result                      = intersection;
                    }
                }
            }
            return result;
        }

        std::optional<Intersection> intersect(const Ray& ray, const RayLimits& limits) noexcept override
        {
            std::optional<Intersection> result;

            auto intermediate_limits{ limits };
            for (const auto& child : m_children) {
                assert(child);
                if (sp::intersect_p(child->m_bounds, ray, intermediate_limits)) {
                    if (const auto intersection = child->intersect(ray, intermediate_limits); intersection) {
                        intermediate_limits.m_t_max = intersection->m_distance;
                        result                      = intersection;
                    }
                }
            }
            return result;
        }

        [[nodiscard]] bool intersect_p(const Ray& ray, const RayLimits& limits) const noexcept override
        {
            for (const auto& child : m_children) {
                assert(child);
                if (sp::intersect_p(child->m_bounds, ray, limits)) {
                    if (child->intersect_p(ray, limits)) {
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

        [[nodiscard]] std::optional<LightIntersection> intersect_lights(const Ray& ray, const RayLimits& limits) noexcept override
        {
            return m_primitives.intersect_lights(ray, limits);
        }

        [[nodiscard]] std::optional<Intersection> intersect(const Ray& ray, const RayLimits& limits) noexcept override
        {
            return m_primitives.intersect(ray, limits);
        }

        [[nodiscard]] bool intersect_p(const Ray& ray, const RayLimits& limits) const noexcept override
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
    [[nodiscard]] std::optional<LightIntersection> intersect_lights_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        assert(m_root);
        return m_root->intersect_lights(ray, limits);
    }

    [[nodiscard]] std::optional<Intersection> intersect_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        assert(m_root);
        return m_root->intersect(ray, limits);
    }

    [[nodiscard]] bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        assert(m_root);
        return m_root->intersect_p(ray, limits);
    }

    [[nodiscard]] BBox3 get_world_bounds_impl() const noexcept override
    {
        return m_root->m_bounds;
    }

    [[nodiscard]] bool is_bounded_impl() const noexcept override
    {
        return true;
    }

    template <typename Iterator>
        requires std::random_access_iterator<Iterator>
    static auto partition(std::size_t dimension, const Point3& absolute_split_location, Iterator first, Iterator last)
    {
        // Put all of the elements that are totally to the left in place
        const auto left =
                std::partition(/*std::execution::par_unseq,*/ first, last, [absolute_split_location, dimension](const auto& a) {
                    return a->get_world_bounds().get_upper()[dimension] < absolute_split_location[dimension];
                });

        // Put all of the elements that are totally to the left in place
        const auto right =
                std::partition(/*std::execution::par_unseq,*/ left, last, [absolute_split_location, dimension](const auto& a) {
                    return a->get_world_bounds().get_lower()[dimension] < absolute_split_location[dimension];
                });
        return std::make_pair(left, right);
    }

    // Reorders elements and returns two additional iterators (called left and right) such that:b
    // Elements in [first, left) are completely to the left of the split location
    // Elements in [right, last) are completely to the right of the split location
    // Elements in [left, right) straddle the split location
    template <typename Iterator>
        requires std::random_access_iterator<Iterator>
    static auto split(std::size_t dimension, const Point3& absolute_split_location, Iterator first, Iterator last)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            const auto [left, right] = partition((dimension + i) % 3, absolute_split_location, first, last);
            if (i == 2 || left != first || right != last) {
                return std::make_pair(left, right);
            }
        }
        return std::make_pair(first, last);
    }

    template <typename Iterator>
        requires std::random_access_iterator<Iterator>
    std::unique_ptr<NodeBase> construct(Iterator first, Iterator last)
    {
        std::unique_ptr<NodeBase> result;
        const auto                num_elements = std::distance(first, last);

        BBox3 bounds;
        std::for_each(std::execution::unseq, first, last, [&bounds](const auto& prim) {
            assert(prim->is_bounded());
            bounds = merge(bounds, prim->get_world_bounds());
        });

        if (num_elements <= k_max_leaf_elements) {
            result.reset(new NodeLeaf{ bounds, first, last });
        } else {
            // TODO: Add surface-area heuristic.
            // Create n (e.g., 16) buckets that span the bounding box in our dimension.
            // For each object, classify which buckets it spans.
            // Calculate SAH for each bucket, and split at the best.
            const auto dimension     = max_dim(bounds.size());
            const auto [left, right] = split(dimension, center(bounds), first, last);

            if (left == first || right == last) {
                result.reset(new NodeLeaf{ bounds, first, last });
                return result;
            }

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
