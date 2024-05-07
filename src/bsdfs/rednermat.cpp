#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class RednerMat final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    RednerMat(const Properties &props) : Base(props) {
        m_albedo = props.texture<Texture>("albedo", .1f);
        m_specular  = props.texture<Texture>("specular", .1f);
        m_roughness = props.texture<Texture>("roughness", .1f);

        m_components.push_back(BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide);
        m_components.push_back(BSDFFlags::GlossyReflection | BSDFFlags::FrontSide);
        m_flags =  m_components[0] | m_components[1];
        dr::set_attr(this, "flags", m_flags);

        // parameters_changed();
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("albedo", m_albedo.get(), +ParamFlags::Differentiable);
        if(m_specular)
            callback->put_object("specular", m_specular.get(), +ParamFlags::Differentiable);
        if(m_roughness)
            callback->put_object("roughness", m_roughness.get(), +ParamFlags::Differentiable);
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float sample1 ,
                                             const Point2f &sample2,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

         bool   has_diffuse  = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0),
                has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);

        Float diffuse_pmf =
             mitsuba::luminance(m_albedo->eval(si, active), si.wavelengths);
         Float specular_pmf =
             mitsuba::luminance(m_specular->eval(si, active), si.wavelengths);
        Float weight_pmf = diffuse_pmf + specular_pmf;
        diffuse_pmf = dr::select(weight_pmf > 0.f, diffuse_pmf / weight_pmf, 0.f);
        specular_pmf = dr::select(weight_pmf > 0.f, specular_pmf / weight_pmf, 0.f);
        // not too sure of this one
        if(dr::neq(has_diffuse, has_specular)){
            diffuse_pmf = dr::select(has_diffuse, 1.0f, 0.f); 
            specular_pmf = dr::select(has_diffuse, 0.0f, 1.f);
        }
        // compute masks Not too sure of those either
        Mask diffuse_mask = diffuse_pmf < sample1;
        Mask specular_mask = !diffuse_mask;

        BSDFSample3f bs   = dr::zeros<BSDFSample3f>();
        bs.eta = 1.0f;

        dr::masked(bs.wo, diffuse_mask) = mitsuba::warp::square_to_cosine_hemisphere(sample2);
        dr::masked(bs.sampled_type, diffuse_mask) = +BSDFFlags::DiffuseReflection;
        dr::masked(bs.sampled_component, diffuse_mask) = 0;

        Float roughness = dr::maximum(m_roughness->eval_1(si, active), 1e-5f);
        Float phong_exponent = roughness_to_phong(roughness);
        // TOOD verify sample2[1]
        Float phi = 2.0 * dr::Pi<Float> * sample2.y();
        Float sin_phi   = dr::sin(phi);
        Float cos_phi = dr::cos(phi);
        Float cos_theta = dr::pow(sample2.x(), dr::rcp(phong_exponent + 2.0f));
        Float sin_theta = dr::safe_sqrt(1.0f - cos_theta * cos_theta);
        Normal3f m =
            Normal3f(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);
        dr::masked(bs.wo, specular_mask) = reflect(si.wi, m);

        dr::masked(bs.sampled_type, specular_mask) = +BSDFFlags::GlossyReflection;
        dr::masked(bs.sampled_component, specular_mask) = 1;

        // Compute PDF
        bs.pdf = pdf(ctx, si, bs.wo, active);
        active &= bs.pdf > 0.f;
        Spectrum result = eval(ctx, si, bs.wo, active);
        return { bs, result / bs.pdf & active };

    }
    Float roughness_to_phong(Float roughness ) const{
        return dr::maximum(2.0f / roughness - 2.0f, 0.0f);
    }
    Float smithG1(Float roughness, Vector3f v) const{
        Float cos_theta = v.z();
        Float tan_theta = dr::safe_sqrt(dr::rcp(cos_theta * cos_theta) - 1.0f);
        Float a         = dr::rcp(dr::sqrt(roughness) * tan_theta);
        Float result = dr::select( a > 1.6, 1.0f, (3.535 * a + 2.181 * a * a) / (1.0f + 2.276 * a + 2.577 * a * a));
        
        return dr::select(dr::eq(tan_theta, 0.0f), 1.0f, result);
    }
    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_diffuse  = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0),
             has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);

        // contribution of the diffuse component
        Vector3f diffuse_contrib =
            dr::select(has_diffuse,
                       m_albedo->eval(si, active) * dr::maximum(wo.z(), 0.0f) /
                           dr::Pi<Float>,
                       0.0f);
        Vector3f m = dr::normalize(si.wi + wo);
        Float roughness = dr::maximum(m_roughness->eval_1(si, active), 1e-5f);
        Float phong_exponent = roughness_to_phong(roughness);
        Vector3f specular_reflectance = m_specular->eval(si, active); // why not eval_1 ?

        Float D = dr::pow(dr::maximum(m.z(), 0.0f), phong_exponent) * (phong_exponent + 2.0f) / (2.0f * dr::Pi<Float>);
        Float G = smithG1(roughness, si.wi) * smithG1(roughness, wo);
        Vector3f F = specular_reflectance + (1.0f - specular_reflectance) * dr::pow(dr::maximum(1.0f - dr::abs( dr::dot(m, wo)), 0.0f), 5.0f);
        Vector3f specular_contrib =
            dr::select(wo.z() > 0.0f, F * D * G / (4.0f * si.wi.z()), 0.0f);
        specular_contrib = dr::select(has_specular, specular_contrib, 0.0f);
        return specular_contrib + diffuse_contrib;
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_diffuse  = ctx.is_enabled(BSDFFlags::DiffuseReflection, 0), 
             has_specular = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);
            

        Float diffuse_pmf =
            mitsuba::luminance(m_albedo->eval(si, active), si.wavelengths);
        Float specular_pmf =
            mitsuba::luminance(m_specular->eval(si, active), si.wavelengths);
        Float weight_pmf = diffuse_pmf + specular_pmf;
        diffuse_pmf = dr::select(weight_pmf > 0.f, diffuse_pmf / weight_pmf, 0.f);
        specular_pmf = dr::select(weight_pmf > 0.f, specular_pmf / weight_pmf, 0.f);
        // not too sure of this one
        if(dr::neq(has_diffuse, has_specular)){
            diffuse_pmf = dr::select(has_diffuse, 1.0f, 0.f); 
            specular_pmf = dr::select(has_diffuse, 0.0f, 1.f);
        }
        // compute diffuse PDF
        Float diffuse_pdf = mitsuba::warp::square_to_cosine_hemisphere_pdf(wo) * diffuse_pmf;
        // compute specular PDF
        Vector3f m = dr::normalize(si.wi + wo);
        Float roughness = dr::maximum(m_roughness->eval_1(si, active), 1e-2f);
        Float phong_exponent = roughness_to_phong(roughness);
        Float D = dr::pow(dr::maximum(m.z(), 0.0f), phong_exponent) * (phong_exponent + 2.0f) / (2.0f * dr::Pi<Float>);
        Float specular_pdf = dr::select(dr::abs(dr::dot(m, wo)) > 0.0f, specular_pmf * D * dr::maximum(m.z(), 0.0f) / (4.0f * dr::abs(dr::dot(m, wo))), 0.0f);
        return dr::select(wo.z() > 0.0f, diffuse_pdf + specular_pdf, 0.0f);
    }

    std::pair<Spectrum, Float> eval_pdf(const BSDFContext &ctx,
                                        const SurfaceInteraction3f &si,
                                        const Vector3f &wo,
                                        Mask active) const override {
        return {eval(ctx, si, wo, active), pdf(ctx, si, wo, active)};
    }

    Spectrum eval_diffuse_reflectance(const SurfaceInteraction3f &si,
                                      Mask active) const override {
        return m_albedo->eval(si, active);
    }
    UInt32 get_texel_index(const SurfaceInteraction3f &si,
            const std::string &reuse_texture,
            Mask active) const override {
       const std::vector<ref<Texture>> textures{ m_albedo, m_specular, m_roughness};

        for (unsigned int i = 0; i < textures.size(); ++i) {
            if (textures[i]->id() == reuse_texture)
                return textures[i]->get_texel_index(si, active);
        }
        return 0u;
    }
    std::string to_string() const override {
        std::ostringstream oss;
        oss << "RednerMat[" << std::endl
            << "  albedo = " << string::indent(m_albedo) << std::endl
            << "  specular = " << string::indent(m_specular) << std::endl
            << "  roughness = " << string::indent(m_roughness) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    ref<Texture> m_albedo;
    ref<Texture> m_specular;
    ref<Texture> m_roughness;
};

MI_IMPLEMENT_CLASS_VARIANT(RednerMat, BSDF)
MI_EXPORT_PLUGIN(RednerMat, "custom bsdf")
NAMESPACE_END(mitsuba)
