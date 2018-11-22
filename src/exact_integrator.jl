struct ExactIntegrator{T,R,V,ET<:AbstractIntegrator{T}} <: AbstractIntegrator{T}
    u::V
    v::V
    E::ET
    function ExactIntegrator(u::V, E::ET) where {T,V,ET<:AbstractIntegrator{T}}
        checkorder(E, ExactIntegrator)
        new{T,real(T),V,ET}(u, similar(u), E)
    end
end

order(::ExactIntegrator) = Inf
order(::Type{ExactIntegrator}) = Inf

function step!(ei::ExactIntegrator{T,R,V,ET}, t::R, Δt::R) where {T,R,V,ET}
    integrate!(ei.E, t, Δt)
    expv!(ei.v, ei.E, Δt, ei.u)
    inplace(ei.E) || copyto!(ei.u, ei.v)
end
