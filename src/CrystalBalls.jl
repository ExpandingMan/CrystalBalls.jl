__precompile__(true)
module CrystalBalls

# abstract types
abstract type AbstractLattice end
abstract type AbstractField end
abstract type AbstractGaugeField <: AbstractField end
export AbstractLattice, AbstractField, AbstractGaugeField

import Base.rand
import Base: getindex

include("unitary.jl")
include("lattice.jl")
include("gaugefields.jl")

end # module
