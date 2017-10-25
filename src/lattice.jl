# here we keep code related to the lattice itself
# note that for the time being I will mostly specialize to 2 dimensions, later I will add code
# for 3 and 4 dimensions. This is being done to avoid splatting.

# TODO should probably implement specialized spacetime vectors using StaticArrays
# TODO fix spacetime indexing scheme

Base.getindex(lat::AbstractLattice, idx::Tuple) = idx
Base.getindex(lat::AbstractLattice, idx...) = getindex(lat, idx)


"""
    fieldsites(::Type{T}, lat::AbstractLattice)

Creates an uninitialized array with elements of type T for every poin on the lattice.
"""
function fieldsites(::Type{T}, lat::AbstractLattice) where T
    Array{T, lat.D+1}((lat.L for i ∈ 1:lat.D)..., lat.T)
end
export fieldsites



"""
    LatticeToroidal <: AbstractLattice

A type for storing parameters of a simple Euclidean rectilinear lattice with D spatial dimensions
and toroidal topology in all directions.
"""
struct LatticeToroidal <: AbstractLattice
    D::Int  # number of spatial dimensions of the lattice
    a::Float64  # lattice spacing in GeV-1
    L::Int  # lattice length in units of a
    T::Int  # lattice duration in units of a
end
export LatticeToroidal

LatticeToroidal(D::Int, a::Float64, L::Int) = LatticeToroidal(D, a, L, L)


nsites(l::AbstractLattice) = l.T*(l.L^l.D)
export nsites


# this is an implementation that doesn't involve splatting
function Base.sub2ind(M::Tuple, m::Tuple)
    1 + sum((m[n]-1)*prod(M[1:(n-1)]) for n ∈ 1:length(M))
end


# helper function for index computers
@inline function _loopedidx(i::Integer, L::Integer)
    # pretty sure this has to be different for negative and positive?
    # Note to self: think in terms of 0-based indexing, otherwise confusing as fuck
    if i < 1
        L + (i % L)
    else
        ((i-1) % L) + 1
    end
end

# these duplications are to avoid splatting
# WARNING right now these don't check lattice dimensions
# for rank-2
function Base.getindex(lat::LatticeToroidal, i::Integer, j::Integer)
    _loopedidx(i, lat.L), _loopedidx(j, lat.T)
end
# for rank-3
function Base.getindex(lat::LatticeToroidal, i::Integer, j::Integer, k::Integer)
    _loopedidx(i, lat.L), _loopedidx(j, lat.L), _loopedidx(k, lat.T)
end
# for rank-4
function Base.getindex(lat::LatticeToroidal, i::Integer, j::Integer, k::Integer, l::Integer)
    _loopedidx(i, lat.L), _loopedidx(j, lat.L), _loopedidx(k, lat.L), _loopedidx(l, lat.L)
end
Base.getindex(lat::LatticeToroidal, idx...) = getindex(lat, idx)
Base.getindex(lat::LatticeToroidal, idx::Tuple) = getindex(lat, idx...)
Base.getindex(lat::LatticeToroidal, idx::SVector) = getindex(lat, Tuple(idx))


# TODO compiler will not know size of these!
spacetimebasisvec(μ::Integer, N::Integer) = SVector{N}([i ≠ μ ? 0 : 1 for i ∈ 1:N])
spacetimebasisvec(μ::Integer, lat::AbstractLattice) = spacetimebasisvec(μ, lat.D + 1)
export spacetimebasisvec

# TODO this is probably expensive and shitty
inc(x::SVector{N}, μ::Integer) where N = SVector{N}([i ≠ μ ? x[i] : x[i]+1 for i ∈ 1:N])
Base.dec(x::SVector{N}, μ::Integer) where N = SVector{N}([i ≠ μ ? x[i] : x[i]-1 for i ∈ 1:N])
export inc
