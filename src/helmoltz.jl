export HelmoltzSolver, solve!, update!

# Helmoltz solver
mutable struct HelmoltzSolver{T, P, Q<:QuasiTridiagonal, V<:Vector{T}}
    Be::Q #
    Bo::Q #
    ge::V #
    go::V #
    l::V  #
    d::V  #
    u::V  #
    function HelmoltzSolver(P::Int, ::Type{T}=Float64) where {T}
        # only available for an even number of expansion coefficients
        isodd(P) || throw(ArgumentError("invalid size P"))

        # this is the size of the problems to be solved
        M = div(P+1, 2)

        # create two banded solver of half the size
        Be = QuasiTridiagonal(M, T)
        Bo = QuasiTridiagonal(M, T)

        # these are the right hand sides
        ge = zeros(T, M)
        go = zeros(T, M)

        _c(p)    = p == 0  ? 2 : 1
        _β(p, P) = p > P-2 ? 0 : 1

        # precompute and store coefficient of the equations from p = 2:P
        l, d, u = zeros(T, P), zeros(T, P), zeros(T, P) 
        for p ∈ 2:P # so that l[p] does not know about julia 1-based indexing
            l[p] = _c(p-2)/(4p*(p-1))
            d[p] = _β(p, P)/2/(p^2 - 1)
            u[p] = _β(p+2, P)/(4*p*(p+1))
        end

        return new{T, P, QuasiTridiagonal{T, M}, Vector{T}}(Be, Bo, ge, go, l, d, u)
    end
end

"""
    Update the banded solvers with new coefficients and factorise.
"""
function update!(h::HelmoltzSolver{T, P}, θ₀::Real, θ₁::Real) where {T, P}
    # size of the two banded systems
    M = div(P+1, 2)

    # update matrix elements
    h.Be.b .= 1
    h.Bo.b .= 1

    # Test case: assume P = 7 and M = (P+1)/2 = 4
    # We have these even and odd coefficients
    # u₀,     u₂,     u₄,     u₆,
    #     u₁,     u₃,     u₅,     u₇
    @simd for i ∈ 1:M-1
        h.Be.l[i] =     -θ₁*h.l[2i]
        h.Bo.l[i] =     -θ₁*h.l[2i+1]
        h.Be.d[i] = θ₀ + θ₁*h.d[2i]
        h.Bo.d[i] = θ₀ + θ₁*h.d[2i+1]
        i < M-1 && (h.Be.u[i] = - θ₁*h.u[2i])
        i < M-1 && (h.Bo.u[i] = - θ₁*h.u[2i+1])
    end

    ul!(h.Bo)
    ul!(h.Be)

    return nothing
end


function solve!(h::HelmoltzSolver{T, P}, f::ChebCoeffs{T, P}, u₊::Real, u₋::Real) where {T, P}
    # size fo the problem
    M = div(P+1, 2)

    # apply boundary conditions
    @inbounds begin
        h.ge[1] = (u₊+u₋)*0.5
        h.go[1] = (u₊-u₋)*0.5

        # we start from 2, after the BC
        @simd for i ∈ 2:M
            p = 2*(i-1) # even modes: write equations for p = 2, 4, 6
            fₚ₊₂ =  p+2 ≥ P ? zero(T) : f[p+2]
            h.ge[i] = h.l[p] * f[p-2] - h.d[p] * f[p] + h.u[p] * fₚ₊₂
        
            p = 2*i-1  # odd modes: write equations for p = 3, 5, 7
            fₚ₊₂ =  p+2 ≥ P ? zero(T) : f[p+2]
            h.go[i] = h.l[p] * f[p-2] - h.d[p] * f[p] + h.u[p] * fₚ₊₂
        end

        # solve two systems
        ldiv!(h.Be, h.ge)
        ldiv!(h.Bo, h.go)

        # copy solution back to input argument f
        @simd for i ∈ 1:M
            f[2i-2] = h.ge[i]
            f[2i-1] = h.go[i]
        end
    end
    return f
end