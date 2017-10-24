# here we keep code related to the lattice itself

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


# holy fuck does this seem inefficient
function Base.getindex(lat::LatticeToroidal, idx::Tuple)
    newidx = 0
    for i ∈ 1:(length(idx)-1)
        if idx[i] > lat.D
            newidx += ((idx[i]-1) % lat.D) + 1
        else
            newidx += idx[i]
        end
    end
    if idx[end] > lat.T
        newidx[end] = ((idx[end]-1) % lat.T) + 1
    else
        newidx[end] = idx[end]
    end
    newidx
end
Base.getindex(lat::LatticeToroidal, idx...) = getindex(lat, idx)

