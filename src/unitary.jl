# code for generating uniform random unitary matrices
# TODO for now we'll just do SU(2)
# TODO generators should probably be of the Julia Hermitian type

# obviously this does not in itself guarantee unitarity
const UnitaryMatrix = Matrix{Complex{Float64}}

# will eventually need to replace with boilerplate code for random SU(N) generators
const σ₁ = Complex{Float64}[0 1; 1 0]
const σ₂ = Complex{Float64}[0 -im; im 0]
const σ₃ = Complex{Float64}[1 0; 0 -1]

const σ⃗ = [σ₁, σ₂, σ₃]

const one2by2 = eye(Complex{Float64}, 2)


"""
    randunitvec(N::Integer)

Creates are random unit vector in Rᴺ.
"""
function randunitvec(::Type{T}, N::Integer) where T
    x = randn(T, N)
    r = vecnorm(x)
    x/r
end
randunitvec(N::Integer) = randunitvec(Float64, N)


"""
    randazimuth()

Creates a random number on [0, 2π).
"""
randazimuth(::Type{T}) where T = 2π*rand(T)
randazimuth() = randazimuth(Float64)


"""
    randSU2_ineff()

An elegant way of generating uniformly distributed random SU(2)
matrices.  This version appears inefficient, as it probably does
lots of dense matrix multiplication.
"""
randSU2_ineff(ϕ::AbstractFloat, n::AbstractVector) = expm(im*ϕ*n'*σ⃗)
randSU2_ineff() = randSU2_ineff(randazimuth(), randunitvec(3))
export randSU2_ineff


"""
    SU2(ϕ::AbstractFloat, n::AbstractVector)

Creat a vector from SU(2) with phase `ϕ*n`.
"""
function SU2(ϕ::AbstractFloat, n::AbstractVector)
    one2by2*cos(ϕ) + im*(n'*σ⃗)*sin(ϕ)
end
export SU2


"""
    randSU2(ϕmax::AbstractFloat)
    randSU2()

Generate a random matrix uniformly distributed on SU(2).
Only involves scalar multiplication and dense matrix addition.

If Φmax is provided, the magnitude of the phase will lie in the
interval [0, ϕmax) (the randomly generated unit vector can be negative).
"""
randSU2(ϕmax::AbstractFloat) = SU2(ϕmax*rand(), randunitvec(3))
randSU2() = randSU2(randazimuth())
export randSU2


"""
    tr(A)

Returns the trace of a matrix. This is simply a convenient alias
for the Julia function `trace`.
"""
tr(A) = trace(A)
export tr
