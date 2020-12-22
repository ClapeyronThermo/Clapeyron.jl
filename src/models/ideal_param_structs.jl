abstract type IdealParams end

struct MonomerParams <: IdealParams
    Mw::Dict
end

struct WalkerParams <: IdealParams
    Mw::Dict
    Nrot::Dict
    theta_V::Dict
    deg_V::Dict
end

struct ReidParams <: IdealParams
    poly_coef::Dict
end

struct WilhoitParams <: IdealParams
    poly_coef::Dict
end

struct NASAParams <: IdealParams
    poly_coef::Dict
end
