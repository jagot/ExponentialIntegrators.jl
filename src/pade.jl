abstract type PadéExponentiator{T} <: AbstractExponentiator{T} end

# This implementation only supports time-independent operators, since
# the implicit step is pre-factorized.
struct CrankNicolson{T,R,F,B,V} <: PadéExponentiator{T}
    fwd::F
    bwd::B
    Δt::R
    tmp::V
end

CrankNicolson(A, Δt::R, μ::Number=one(R)) where R =
    CrankNicolson(I + μ*Δt*A/2, factorize(I - μ*Δt*A/2),
                  Vector{typeof(μ)}(undef, size(A,1)))

# Since the implicit step is pre-factorized, it is not possible to
# integrate over the time-step. Returning `nothing` signifies this.
operator(::CrankNicolson) = nothing

inplace(::CrankNicolson) = true
order(::CrankNicolson) = 2

function expv!(_, CN::CrankNicolson{T,R}, Δt::R, u)
    @assert Δt == CN.Δt
    mul!(CN.tmp, CN.fwd, u)
    ldiv!(u, CN.bwd, CN.tmp)
end
