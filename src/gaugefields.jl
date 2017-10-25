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
struct GaugeFieldSU2{n} <: AbstractGaugeField{n}
    field::Array{WilsonLineVar, n}
    lattice::AbstractLattice

    GaugeFieldSU2{n}(st::AbstractArray{WilsonLineVar,n}, lat::AbstractLattice) where n = new(st, lat)
end
export GaugeFieldSU2

function GaugeFieldSU2(spacetime::AbstractArray{WilsonLineVar,n}, lat::AbstractLattice) where n
    if lat.D + 1 ≠ n
        ArgumentError("Incompatible spacetime dimensions.")
    end
    GaugeFieldSU2{n}(spacetime, lat)
end

# for now this just creates a "vacuum" with A=0 everywhere
function GaugeFieldSU2(lat::AbstractLattice)
    # note inefficiency due to splatting
    spacetime = fieldsites(WilsonLineVar, lat)
    for idx ∈ eachindex(spacetime)
        spacetime[idx] = wilsonlineunity(2, lat)
    end
    GaugeFieldSU2(spacetime, lat)
end


# TODO ensure that splatting really is ok here!
# gets the gauge field at spacetime location x1, t
function Base.getindex(f::GaugeFieldSU2, x::Tuple)
    idx = getindex(f.lattice, x)
    getindex(f.field, idx...)
end
Base.getindex(f::GaugeFieldSU2, x::SVector) = getindex(f, Tuple(x))
Base.getindex(f::GaugeFieldSU2, x...) = getindex(f, x...)


# TODO should do this for abstract gauge field, must make parametric
# this gets each spatial index of the gauge field
function eachspacetimeindex(U::AbstractGaugeField{n}) where n
    (SVector{n}(ind2sub(size(U.field), x)) for x ∈ eachindex(U.field))
end
export eachspacetimeindex


function Base.rand(::Type{GaugeFieldSU2}, lat::AbstractLattice)
    spacetime = fieldsites(WilsonLineVar, lat)
    for idx ∈ eachindex(spacetime)
        spacetime[idx] = randwilsonlineSU2(lat)
    end
    GaugeFieldSU2(spacetime, lat)
end


"""
    untracedplaquette(μ::Int, ν::Int, x::Tuple)

Compute the value of a plaquetee before taking the trace.  In other words
\$\$ U_{\mu}(x)U_{\nu}(x+e_{\mu})U^{\dagger}_{\mu}(x+e_{\nu})U^{\dagger}_{\nu}(x) \$\$
"""
function untracedplaquette(U::AbstractGaugeField, μ::Integer, ν::Integer, x::SVector)
    U[x][μ]*U[inc(x,μ)][ν]*U[inc(x,ν)][μ]'*U[x][ν]'
end
export untracedplaquette

"""
    plaquette(U::AbstractGaugeField, μ::Integer, ν::Integer, x::SVector)

Compute the value of the plaquette at location `x` over the gauge field.
"""
plaquette(U::AbstractGaugeField, μ::Integer, ν::Integer, x::SVector) = tr(untracedplaquette(U,μ,ν,x))
export plaquette


"""
    wilsonlagrangian(U::AbstractGaugeField, β::AbstractFloat, x::SVector)

Compute the lowest order Wilson plaquette Lagrangian at location `x` and inverse temperature `β`.

Note that this has a factor of a^4 in front of it relative to the usual G^2 Yang-Mills kinetic term.
"""
function wilsonlagrangian(U::AbstractGaugeField, β::AbstractFloat, x::SVector)
    ℒ = 0.0
    for μ ∈ 1:(U.lattice.D+1)
        for ν ∈ 1:(μ-1)
            ℒ += 1.0 - real(plaquette(U, μ, ν, x))/3.0
        end
    end
    β*ℒ
end
export wilsonlagrangian


"""
    wilsonaction(U::AbstractGaugeField, β::AbstractFloat)

Compute the Wilson plaquette action at inverse temperature `β`.
"""
function wilsonaction(U::AbstractGaugeField, β::AbstractFloat)
    S = 0.0
    for x ∈ eachspacetimeindex(U)
        S += wilsonlagrangian(U, β, x)
    end
    S/a^(3 - U.lattice.D) # this factor to fix the powers of a
end
export wilsonaction
