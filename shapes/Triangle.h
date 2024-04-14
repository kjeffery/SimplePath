#pragma once

/// @author Keith Jeffery

#include "Shape.h"

#include "../math/Vector3.h"

#include <algorithm>
#include <array>
#include <execution>
#include <memory>
#include <utility>
#include <vector>

namespace sp {
class Mesh;

void log_extents(const Mesh&, const char* const str);

class Mesh
{
public:
    Mesh(std::vector<std::size_t> indices,
         std::vector<Point3>      vertices,
         std::vector<Normal3>     normals,
         const AffineSpace&       object_to_world)
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
    : Shape(AffineSpace::identity(), AffineSpace::identity())
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
        const Point3& p0 = m_mesh->m_vertices[m_indices[0]];
        const Point3& p1 = m_mesh->m_vertices[m_indices[1]];
        const Point3& p2 = m_mesh->m_vertices[m_indices[2]];

        const float A = p0.x - p1.x;
        const float B = p0.y - p1.y;
        const float C = p0.z - p1.z;

        const float D = p0.x - p2.x;
        const float E = p0.y - p2.y;
        const float F = p0.z - p2.z;

        const float G = ray.get_direction().x;
        const float H = ray.get_direction().y;
        const float I = ray.get_direction().z;

        const float J = p0.x - ray.get_origin().x;
        const float K = p0.y - ray.get_origin().y;
        const float L = p0.z - ray.get_origin().z;

        const float EIHF = msub(E, I, H * F); // E*I - H*F
        const float GFDI = msub(G, F, D * I); // G*F - D*I
        const float DHEG = msub(D, H, E * G); // D*H - E*G

        const float denom = madd(A, EIHF, madd(B, GFDI, C * DHEG)); // A*EIHF + B*GFDI + C*DHEG
        if (denom == 0) {
            return {};
        }

        const float beta = madd(J, EIHF, madd(K, GFDI, L * DHEG)) / denom; // J*EIHF + K*GFDI + L*DHEG
        if (beta <= 0.0f || beta >= 1.0f) {
            return {};
        }

        const float AKJB = msub(A, K, J * B); // A*K - J*B
        const float JCAL = msub(J, C, A * L); // J*C - A*L
        const float BLKC = msub(B, L, K * C); // B*L - K*C

        const float gamma = madd(I, AKJB, madd(H, JCAL, G * BLKC)) / denom; // I*AKJB + H*JCAL + G*BLKC
        if (gamma <= 0.0f || beta + gamma >= 1.0f) {
            return {};
        }

        const float t = -madd(F, AKJB, madd(E, JCAL, D * BLKC)) / denom; // F*AKJB + E*JCAL + D*BLKC
        if (t < limits.m_t_min || t > limits.m_t_max) {
            return {};
        }

        const Normal3& n0 = m_mesh->m_normals[m_indices[0]];
        const Normal3& n1 = m_mesh->m_normals[m_indices[1]];
        const Normal3& n2 = m_mesh->m_normals[m_indices[2]];

        // Use the barycentric coordinates to interpolate our normal.
        const float alpha = 1.0f - beta - gamma;
        // const Normal3 normal = normalize(alpha * n0 + beta * n1 + gamma * n2);
        const Normal3 normal = normalize(madd(alpha, n0, madd(beta, n1, gamma * n2)));

        Intersection isect;
        isect.m_point    = ray(t);
        isect.m_normal   = normal;
        isect.m_distance = t;

        return { isect };
    }

    [[nodiscard]] bool intersect_p_impl(const Ray& ray, const RayLimits& limits) const noexcept override
    {
        const Point3& p0 = m_mesh->m_vertices[m_indices[0]];
        const Point3& p1 = m_mesh->m_vertices[m_indices[1]];
        const Point3& p2 = m_mesh->m_vertices[m_indices[2]];

        const float A = p0.x - p1.x;
        const float B = p0.y - p1.y;
        const float C = p0.z - p1.z;

        const float D = p0.x - p2.x;
        const float E = p0.y - p2.y;
        const float F = p0.z - p2.z;

        const float G = ray.get_direction().x;
        const float H = ray.get_direction().y;
        const float I = ray.get_direction().z;

        const float J = p0.x - ray.get_origin().x;
        const float K = p0.y - ray.get_origin().y;
        const float L = p0.z - ray.get_origin().z;

        const float EIHF = msub(E, I, H * F); // E*I - H*F
        const float GFDI = msub(G, F, D * I); // G*F - D*I
        const float DHEG = msub(D, H, E * G); // D*H - E*G

        const float denom = madd(A, EIHF, madd(B, GFDI, C * DHEG)); // A*EIHF + B*GFDI + C*DHEG
        if (denom == 0) {
            return false;
        }

        const float beta = madd(J, EIHF, madd(K, GFDI, L * DHEG)) / denom; // J*EIHF + K*GFDI + L*DHEG
        if (beta <= 0.0f || beta >= 1.0f) {
            return false;
        }

        const float AKJB = msub(A, K, J * B); // A*K - J*B
        const float JCAL = msub(J, C, A * L); // J*C - A*L
        const float BLKC = msub(B, L, K * C); // B*L - K*C

        const float gamma = madd(I, AKJB, madd(H, JCAL, G * BLKC)) / denom; // I*AKJB + H*JCAL + G*BLKC
        if (gamma <= 0.0f || beta + gamma >= 1.0f) {
            return false;
        }

        const float t = -madd(F, AKJB, madd(E, JCAL, D * BLKC)) / denom; // F*AKJB + E*JCAL + D*BLKC
        if (t < limits.m_t_min || t > limits.m_t_max) {
            return false;
        }

        return true;
    }

    [[nodiscard]] BBox3 get_object_bounds() const noexcept override
    {
        // Vertices are transformed into world bounds at mesh creation.
        BBox3 bounds;
        for (std::size_t i = 0; i < 3; ++i) {
            const Point3& p = m_mesh->m_vertices[m_indices[i]];
            bounds.extend(get_world_to_object()(p));
        }
        return bounds;
    }

    [[nodiscard]] BBox3 get_world_bounds_impl() const noexcept override
    {
        // Vertices are already transformed into world bounds at mesh creation.

        BBox3 bounds;
        for (std::size_t i = 0; i < 3; ++i) {
            const Point3& p = m_mesh->m_vertices[m_indices[i]];
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
