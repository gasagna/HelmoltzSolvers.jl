using StaticArrays

export CoupledHelmoltzSolver
# Solve the problem
#
#       / θ₀ u''(y) - θ₁ u(y)        = r(y)
#      |  θ₂ v''(y) - θ₃ v(y) - u(y) = 0
#       \ v(±1) = v'(±1) = 0
#
# using the Chebychev tau method with an influence matrix technique,
# and return the solution `v`, overwriting the input argument `r`.

struct CoupledHelmoltzSolver{T, P, H<:HelmoltzSolver{T, P}, C<:ChebCoeffs{T, P}}
    hu::H
    hv::H
    vₛ::NTuple{3, C}
    function CoupledHelmoltzSolver(P::Int, ::Type{T}=Float64) where {T}
        hu = HelmoltzSolver(P, T)
        hv = HelmoltzSolver(P, T)
        vₛ = ntuple(i->ChebCoeffs(P, T), 3)
        return new{T, P, typeof(hu), typeof(vₛ[1])}(hu, hv, vₛ)
    end
end

function solve!(solver::CoupledHelmoltzSolver{T, P},
                    θs::NTuple{4, Real},
                     r::ChebCoeffs{T, P}) where {T, P}

    # expand coeffs
    θ₀, θ₁, θ₂, θ₃ = θs

    # aliases for the partial solutions
    vₚ, v₊, v₋ = solver.vₛ

    # update quasi-tridiagonal solvers with coefficients
    # this also trigger the UL factorisation of the systems
    update!(solver.hu, θ₀, θ₁)
    update!(solver.hv, θ₂, θ₃)

    # `vₚ` and will be overwritten with the solution of the inhomogeneous
    # B problem and eventually with the full solution `v`.
    vₚ.data .= r.data
    v₊.data .= 0
    v₋.data .= 0

    # ~~~~ Solve the three different problems ~~~
    # inhomogeneous B problem
    solve!(solver.hu, vₚ, 0, 0)
    solve!(solver.hv, vₚ, 0, 0)

    # homogeneous B+ problem
    solve!(solver.hu, v₊, 1, 0)
    solve!(solver.hv, v₊, 0, 0)

    # homogeneous b_ problem
    solve!(solver.hu, v₋, 0, 1)
    solve!(solver.hv, v₋, 0, 0)

    # ~~~~ Influence matrix equations ~~~
    #  Note that julia is column major
    A = SMatrix{2, 2}( _ddy(v₊, Val(:right)),  _ddy(v₊, Val(:left)), _ddy(v₋, Val(:right)), _ddy(v₋, Val(:left)))
    b = SVector{2}(-_ddy(vₚ, Val(:right)), -_ddy(vₚ, Val(:left)))
    δ₊, δ₋ = A\b

    # construct full solution
    r.data .= vₚ.data .+ δ₊ .* v₊.data .+ δ₋ .* v₋.data

    return r
end