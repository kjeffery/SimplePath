#pragma once

/// @author Keith Jeffery

#include "Shape.h"

#include "../math/Vector3.h"

#include <algorithm>
#include <execution>
#include <memory>
#include <utility>
#include <vector>

#include "../math/Transformation.h"

namespace sp {
class Mesh;

void log_extents(const Mesh&, const char* str);

class Mesh
{
public:
    Mesh(std::vector<std::size_t>    indices,
         std::vector<Point3>         vertices,
         std::vector<Normal3>        normals,
         const AffineTransformation& object_to_world)
    : m_indices(std::move(indices))
    , m_vertices(std::move(vertices))
    , m_normals(std::move(normals))
    {
        log_extents(*this, "Pre-transform");

        // Unlike other objects, we don't transform the ray at intersection test time: we pre-transform all of our
        // triangles.
        std::transform( /*std::execution::par_unseq,*/
            m_vertices.cbegin(),
            m_vertices.cend(),
            m_vertices.begin(),
            [&object_to_world](auto& v) { return object_to_world(v); });

        std::transform( /*std::execution::par_unseq,*/
            m_normals.cbegin(),
            m_normals.cend(),
            m_normals.begin(),
            [&object_to_world](auto& n) { return object_to_world(n); });

        log_extents(*this, "Post-transform");
    }

    Mesh(Mesh&&)                 = default;
    Mesh(const Mesh&)            = default;
    Mesh& operator=(Mesh&&)      = default;
    Mesh& operator=(const Mesh&) = default;

    [[nodiscard]] std::size_t get_num_triangles() const noexcept
    {
        const std::size_t s = m_indices.size();
        assert(s % 3u == 0u);
        return s / 3u;
    }

    std::vector<std::size_t> m_indices;
    std::vector<Point3>      m_vertices;
    std::vector<Normal3>     m_normals;
};

class Triangle : public Shape
{
public:
#if 0
    Triangle(AffineSpace           object_to_world,
             AffineSpace           world_to_object,
             std::shared_ptr<Mesh> mesh,
             std::size_t           triangle_number)
    : Shape(object_to_world, world_to_object)
    , m_mesh{ mesh }
    , m_indices{ mesh->m_indices.data() + triangle_number * 3 }
    {
    }
#else
    Triangle(std::shared_ptr<Mesh> mesh, std::size_t triangle_number)
    : Shape()
    , m_mesh{ std::move(mesh) }
    , m_indices{ m_mesh->m_indices.data() + triangle_number * 3 }
    {
    }
#endif

private:
    [[nodiscard]] std::optional<LightIntersection> intersect_lights_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        return {};
    }

    [[nodiscard]] std::optional<Intersection> intersect_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        const auto& p0 = m_mesh->m_vertices[m_indices[0]];
        const auto& p1 = m_mesh->m_vertices[m_indices[1]];
        const auto& p2 = m_mesh->m_vertices[m_indices[2]];

        const auto A = p0.x - p1.x;
        const auto B = p0.y - p1.y;
        const auto C = p0.z - p1.z;

        const auto D = p0.x - p2.x;
        const auto E = p0.y - p2.y;
        const auto F = p0.z - p2.z;

        const auto G = ray.get_direction().x;
        const auto H = ray.get_direction().y;
        const auto I = ray.get_direction().z;

        const auto J = p0.x - ray.get_origin().x;
        const auto K = p0.y - ray.get_origin().y;
        const auto L = p0.z - ray.get_origin().z;

        const auto EIHF = msub(E, I, H * F); // E*I - H*F
        const auto GFDI = msub(G, F, D * I); // G*F - D*I
        const auto DHEG = msub(D, H, E * G); // D*H - E*G

        const auto denom = madd(A, EIHF, madd(B, GFDI, C * DHEG)); // A*EIHF + B*GFDI + C*DHEG
        if (denom == 0) {
            return {};
        }

        const auto beta = madd(J, EIHF, madd(K, GFDI, L * DHEG)) / denom; // J*EIHF + K*GFDI + L*DHEG
        if (beta <= 0.0f || beta >= 1.0f) {
            return {};
        }

        const auto AKJB = msub(A, K, J * B); // A*K - J*B
        const auto JCAL = msub(J, C, A * L); // J*C - A*L
        const auto BLKC = msub(B, L, K * C); // B*L - K*C

        const auto gamma = madd(I, AKJB, madd(H, JCAL, G * BLKC)) / denom; // I*AKJB + H*JCAL + G*BLKC
        if (gamma <= 0.0f || beta + gamma >= 1.0f) {
            return {};
        }

        const auto t = -madd(F, AKJB, madd(E, JCAL, D * BLKC)) / denom; // F*AKJB + E*JCAL + D*BLKC
        if (t < limits.m_t_min || t > limits.m_t_max) {
            return {};
        }

        const auto& n0 = m_mesh->m_normals[m_indices[0]];
        const auto& n1 = m_mesh->m_normals[m_indices[1]];
        const auto& n2 = m_mesh->m_normals[m_indices[2]];

        // Use the barycentric coordinates to interpolate our normal.
        const auto alpha = 1.0f - beta - gamma;
        // const Normal3 normal = normalize(alpha * n0 + beta * n1 + gamma * n2);
        const auto normal = normalize(madd(alpha, n0, madd(beta, n1, gamma * n2)));

        Intersection isect;
        isect.m_point    = ray(t);
        isect.m_normal   = normal;
        isect.m_distance = t;

        return { isect };
    }

    [[nodiscard]] bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        const auto& p0 = m_mesh->m_vertices[m_indices[0]];
        const auto& p1 = m_mesh->m_vertices[m_indices[1]];
        const auto& p2 = m_mesh->m_vertices[m_indices[2]];

        const auto A = p0.x - p1.x;
        const auto B = p0.y - p1.y;
        const auto C = p0.z - p1.z;

        const auto D = p0.x - p2.x;
        const auto E = p0.y - p2.y;
        const auto F = p0.z - p2.z;

        const auto G = ray.get_direction().x;
        const auto H = ray.get_direction().y;
        const auto I = ray.get_direction().z;

        const auto J = p0.x - ray.get_origin().x;
        const auto K = p0.y - ray.get_origin().y;
        const auto L = p0.z - ray.get_origin().z;

        const auto EIHF = msub(E, I, H * F); // E*I - H*F
        const auto GFDI = msub(G, F, D * I); // G*F - D*I
        const auto DHEG = msub(D, H, E * G); // D*H - E*G

        const auto denom = madd(A, EIHF, madd(B, GFDI, C * DHEG)); // A*EIHF + B*GFDI + C*DHEG
        if (denom == 0) {
            return false;
        }

        const auto beta = madd(J, EIHF, madd(K, GFDI, L * DHEG)) / denom; // J*EIHF + K*GFDI + L*DHEG
        if (beta <= 0.0f || beta >= 1.0f) {
            return false;
        }

        const auto AKJB = msub(A, K, J * B); // A*K - J*B
        const auto JCAL = msub(J, C, A * L); // J*C - A*L
        const auto BLKC = msub(B, L, K * C); // B*L - K*C

        const auto gamma = madd(I, AKJB, madd(H, JCAL, G * BLKC)) / denom; // I*AKJB + H*JCAL + G*BLKC
        if (gamma <= 0.0f || beta + gamma >= 1.0f) {
            return false;
        }

        const auto t = -madd(F, AKJB, madd(E, JCAL, D * BLKC)) / denom; // F*AKJB + E*JCAL + D*BLKC
        if (t < limits.m_t_min || t > limits.m_t_max) {
            return false;
        }

        return true;
    }

    [[nodiscard]] BBox3 get_object_bounds_impl() const noexcept override
    {
        // Vertices are transformed into world bounds at mesh creation.
        BBox3 bounds;
        for (std::size_t i = 0; i < 3; ++i) {
            const auto& p = m_mesh->m_vertices[m_indices[i]];
            bounds.extend(p);
        }
        return bounds;
    }

    [[nodiscard]] BBox3 get_world_bounds_impl() const noexcept override
    {
        // Vertices are already transformed into world bounds at mesh creation.
        BBox3 bounds;
        for (std::size_t i = 0; i < 3; ++i) {
            const auto& p = m_mesh->m_vertices[m_indices[i]];
            bounds.extend(p);
        }
        return bounds;
    }

    [[nodiscard]] bool is_bounded_impl() const noexcept override
    {
        return true;
    }

    std::shared_ptr<Mesh> m_mesh;
    std::size_t*          m_indices;
};
} // namespace sp
