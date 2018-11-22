# By default, an AbstractOperator is assumed to be time-independent.
midpoint!(::Union{AbstractOperator,Nothing}, t, Δt) = nothing
integrate!(::Union{AbstractOperator,Nothing}, t, Δt) = nothing
