# code for gauge fields

# TODO for now the wilson lines will be `Vector`s of SU(2) matrices
# TODO looks like for arbitrary dimensions we'll have lots of splatting inefficiency...
#   this should be ok because most of this stuff should just be for initialization...
#   in the future can fix splatting inefficiency with macros

const WilsonLineVar = Vector{UnitaryMatrix}
export WilsonLineVar

wilsonlineunity(N::Integer, D::Integer) = UnitaryMatrix[eye(Complex{Float64}, N) for d ∈ 1:(D+1)]
wilsonlineunity(N::Integer, lat::AbstractLattice) = wilsonlineunity(N, lat.D)
export wilsonlineunity

randwilsonlineSU2(D::Integer) = UnitaryMatrix[randSU2() for d ∈ 1:(D+1)]
randwilsonlineSU2(lat::AbstractLattice) = randwilsonlineSU2(lat.D)
export randwilsonlineSU2

"""
    GaugeFieldSU2{n} <: AbstractGaugeField

An SU(2) gauge field.  Note that we take the fundamental field variables to be Wilson lines emanating
from lattice site x in direction μ. Here `n` is the total number of dimensions `D+1`.
"""
struct GaugeFieldSU2{n} <: AbstractGaugeField
    field::Array{WilsonLineVar, n}

    GaugeFieldSU2{n}(st::AbstractArray{WilsonLineVar,n}) where n = new(st)
end
export GaugeFieldSU2

GaugeFieldSU2(spacetime::AbstractArray{WilsonLineVar,n}) where n = GaugeFieldSU2{n}(spacetime)

# for now this just creates a "vacuum" with A=0 everywhere
function GaugeFieldSU2(lat::AbstractLattice)
    # note inefficiency due to splatting
    spacetime = fieldsites(WilsonLineVar, lat)
    for idx ∈ eachindex(spacetime)
        spacetime[idx] = wilsonlineunity(2, lat)
    end
    GaugeFieldSU2(spacetime)
end

Base.getindex(f::GaugeFieldSU2, idx...) = getindex(f.field, idx...)


function Base.rand(::Type{GaugeFieldSU2}, lat::AbstractLattice)
    spacetime = fieldsites(WilsonLineVar, lat)
    for idx ∈ eachindex(spacetime)
        spacetime[idx] = randwilsonlineSU2(lat)
    end
    GaugeFieldSU2(spacetime)
end


"""
    untracedplaquette(μ::Int, ν::Int, x::Tuple)

Compute the value of a plaquetee before taking the trace.  In other words
$$ U_{\mu}(x)U_{\nu}(x+e_{\mu})U^{\dagger}_{\mu}(x+e_{\nu})U^{\dagger}_{\nu}(x) $$
"""
function untracedplaquette(field::AbstractGaugeField, μ::Int, ν::Int, x::Tuple)

end
