abstract type KrylovExponentiator{T} <: AbstractExponentiator{T} end

# This is of course only true if the Krylov subspace is large enough.
exact(::KrylovExponentiator) = true
order(::KrylovExponentiator) = Inf

mutable struct Lanczos{T,U,O<:AbstractOperator{T},scale} <: KrylovExponentiator{T}
    A::O
    μ::T
    Ks::KrylovSubspace{U,T,U}
    cache::StegrCache{T,U}
    atol::U
    rtol::U
end

function LanczosMap(::Type{T}, A::O, μ::T, m::Integer; scale::Bool=true, atol=1e-8, rtol=√atol) where {T,O<:AbstractOperator}
    n = size(A,2)
    Ks = KrylovSubspace{T,real(T)}(n,m)
    LanczosMap{T,real(T),O,scale}(A, μ, Ks, StegrCache(T, m), atol, rtol)
end

LanczosMap(::Type{T}, A, m::Integer; kwargs...) where T =
    LanczosMap(T, A, one(T), m; kwargs...)

operator(L::LanczosMap) = L.A

# Scale the operator by the time-step, to reduce the spectral radius.
expv!(v, L::LanczosMap{T,U,true}, Δt::U, u, _) where {T,U} =
    expv!(v, L.μ, Δt*L.A, u, L.Ks, L.cache, atol=L.atol, rtol=L.rtol)

# Apply the time-step when calculating the subspace exponential.
expv!(v, L::LanczosMap{T,U,false}, Δt::U, u, _) where {T,U} =
    expv!(v, L.μ*Δt, L.A, u, L.Ks, L.cache, atol=L.atol, rtol=L.rtol)
