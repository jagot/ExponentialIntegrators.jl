module ExponentialIntegrators

using ExponentialUtilities

import ExponentialUtilities: expv!

abstract type AbstractIntegrator{T} end

order(::AbstractIntegrator) = 1
order(::Type{AbstractIntegrator}) = 1

abstract type AbstractExponentiator{T} end

exact(::AbstractExponentiator) = false
inplace(::AbstractExponentiator) = false
order(::AbstractExponentiator) = 1

# The last argument to `expv!`, the integer, signifies if the
# exponentiation should be performed in the forward direction (1), or
# backward direction (0). This is useful if `E` split further, and is
# part of a symmetric splitting scheme. If the argument is 0, the
# order does not matter (as is the case for the middle term of a
# certain splitting).
expv!(v, E::AbstractExponentiator, Δt, u) = expv!(v, E, Δt, u, 0)

abstract type AbstractOperator{T} end
# Time-dependent operators could be of the kind A(t), A + f(t)B, etc,
# for which the different temporal integration routines would have to
# be special-cased.

checkorder(o,integrator) =
    o < order(integrator) && @warn("Order of $(typeof(o)) lower than required for $(integrator)")

macro swap!(x,y)
    quote
        local tmp = $(esc(x))
        $(esc(x)) = $(esc(y))
        $(esc(y)) = tmp
    end
end

# Temporal integration of operators
include("integrate.jl")

# Exponentiators for single operators
include("krylov.jl")
include("pade.jl")

# Integrators composed of multiple exponentiators
include("strang_splitting.jl")
include("exact_integrator.jl")

export step!, StrangSplitting, ExactIntegrator

end # module
