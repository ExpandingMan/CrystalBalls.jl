__precompile__(true)
module CrystalBalls

using StaticArrays

# abstract types
abstract type AbstractLattice end
# the type parameters give the spacetime dimensionality
abstract type AbstractField{n} end
abstract type AbstractGaugeField{n} <: AbstractField{n} end
export AbstractLattice, AbstractField, AbstractGaugeField

import Base.rand
import Base: getindex, sub2ind, dec

include("unitary.jl")
include("lattice.jl")
include("gaugefields.jl")
include("metropolis.jl")

end # module
