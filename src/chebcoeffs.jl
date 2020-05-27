export ChebCoeffs

# Chebychev series coefficients
struct ChebCoeffs{T, P, V<:AbstractVector{T}} <: AbstractVector{T}
    data::V

    # construct from size and type
    ChebCoeffs(P::Int, ::Type{T}=Float64) where {T} =
        new{T, P, Vector{T}}(zeros(P+1))

    # consttruct from data directly
    ChebCoeffs(v::V) where {T, V<:AbstractVector{T}} =
        new{T, length(v)-1, V}(v)
end


@inline Base.LinearIndices(s::ChebCoeffs) = axes(s, 1)
@inline Base.axes(s::ChebCoeffs{T, P}, d) where {T, P} =
    d == 1 ? (0:P) : Base.OneTo(1)
@inline Base.size(s::ChebCoeffs{T, P}) where {T, P} = (P+1, )

Base.similar(s::ChebCoeffs) = ChebCoeffs(similar(s.data))
Base.copy(s::ChebCoeffs) = ChebCoeffs(copy(s.data))
Base.parent(s::ChebCoeffs) = s.data

Base.@propagate_inbounds function Base.getindex(s::ChebCoeffs, i::Int)
    @boundscheck checkbounds(s.data, i+1)
    @inbounds ret = s.data[i+1]
    return ret
end

Base.@propagate_inbounds function Base.setindex!(s::ChebCoeffs, v, i::Int)
    @boundscheck checkbounds(s.data, i+1)
    @inbounds s.data[i+1] = v
    return v
end

# For fast broadcasting: ref https://discourse.julialang.org/t/why-is-there-a-performance-hit-on-broadcasting-with-offsetarrays/32194
Base.dataids(s::ChebCoeffs) = Base.dataids(parent(s))
Broadcast.broadcast_unalias(dest::ChebCoeffs, src::ChebCoeffs) = 
    parent(dest) === parent(src) ? src : Broadcast.unalias(dest, src)

# Derivative at S = ± 1
function _ddy(s::ChebCoeffs{T, P}, ::Val{:left}) where {T, P}
    @inbounds begin
        dsdy = -s[1]
        @simd for p ∈ 2:2:P
            dsdy +=     p^2 * s[p]
            dsdy -= (p+1)^2 * s[p+1]
        end
    end
    return -dsdy
end

function _ddy(s::ChebCoeffs{T, P}, ::Val{:right}) where {T, P}
    @inbounds begin
        dsdy = s[1]
        @simd for p ∈ 2:P
            dsdy += p^2 * s[p]
        end
    end
    return dsdy
end