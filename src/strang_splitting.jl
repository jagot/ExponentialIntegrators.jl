struct StrangSplitting{T,R,V} <: AbstractIntegrator{T}
    uv::Tuple{V,V}
    E::Vector{<:AbstractExponentiator{T}}
    function StrangSplitting(u::V, E::AbstractExponentiator{T}...) where {T,V}
        foreach(e -> checkorder(e,StrangSplitting), E)
        new{T,real(T),V}((u,similar(u)), Vector{AbstractExponentiator{T}}(E...))
    end
end

order(::StrangSplitting) = 2
order(::Type{StrangSplitting}) = 2

function step!(strang::StrangSplitting{T,R,V}, t::R, Δt::R) where {T,R,V}
    # Integrate all operators over time-step, to second order since
    # the Strang splitting itself introduces a second-order splitting
    # error.
    foreach(e -> midpoint!(operator(e), t, Δt), strang.E)

    a,b = 1,2

    # When the operator is L = A + B + C + D + ..., the Strang
    # splitting is ... (D/2) (C/2) (B/2) A (B/2) (C/2) (D/2) ...

    # ... (D/2) (C/2) (B/2)
    for E in reverse(strang.E[2:end])
        expv!(strang.uv[b], E, Δt/2, strang.uv[a], -1)
        inplace(E) || @swap! a b
    end

    # A
    expv!(strang.uv[b], E₁, Δt, strang.uv[a], 0)
    inplace(strang.E[1]) || @swap! a b

    # (B/2) (C/2) (D/2) ...
    for E in strang.E[2:end]
        expv!(strang.uv[b], E, Δt/2, strang.uv[a], 1)
        inplace(E) || @swap! a b
    end

    # If the propagated solution is in the v vector, copy it to the u
    # vector.
    a == 2 && copyto!(strang.uv[1], strang.uv[2])
end
