#pragma once

#include "Distribution1D.h"
#include "Point2i.h"
#include "../Image/Image.h"

class Distribution2D
{
public:
    explicit Distribution2D(const sp::Image& image)
    : p_conditional{ create_conditional(image) }
    , p_marginal{ create_marginal(p_conditional) }
    {
    }

    [[nodiscard]] auto sample_continuous(const sp::Point2& u) const -> std::pair<sp::Point2, float>
    {
        std::array<float, 2> pdfs;
        std::size_t          v;
        const auto           d1 = p_marginal.sample_continuous(u[1], &pdfs[1], &v);
        const auto           d0 = p_conditional[v].sample_continuous(u[0], &pdfs[0]);
        return { { d0, d1 }, pdfs[0] * pdfs[1] };
    }

    [[nodiscard]] auto pdf(const sp::Point2& p) const -> float
    {
        const auto iu = std::clamp(static_cast<std::size_t>(p[0] * static_cast<float>(p_conditional[0].size())), std::size_t{ 0 },
                                   p_conditional[0].size() - std::size_t{ 1 });
        const auto iv = std::clamp(static_cast<std::size_t>(p[1] * static_cast<float>(p_marginal.size())), std::size_t{ 0 },
                                   p_marginal.size() - std::size_t{ 1 });
        return p_conditional[iv].eval(iu) / p_marginal.integral();
    }

private:
    [[nodiscard]] static auto create_conditional(const sp::Image& image) -> std::vector<Distribution1D>
    {
        std::vector<Distribution1D> conditional;
        conditional.reserve(image.height());
        for (std::size_t v = 0; v < image.height(); ++v) {
            std::vector<float> values;
            values.reserve(image.width());
            for (std::size_t u = 0; u < image.width(); ++u) {
                values.push_back(sp::relative_luminance(image(u, v)));
            }
            conditional.emplace_back(std::move(values));
        }
        return conditional;
    }

    [[nodiscard]] static auto create_marginal(const std::vector<Distribution1D>& conditionals) -> Distribution1D
    {
        std::vector<float> marginal_values;
        for (const auto& conditional : conditionals) {
            marginal_values.push_back(conditional.integral());
        }
        return Distribution1D{ std::move(marginal_values) };
    }

    std::vector<Distribution1D> p_conditional;
    Distribution1D              p_marginal;
};
