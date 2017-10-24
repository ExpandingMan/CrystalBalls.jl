# code for generating uniform random unitary matrices
# TODO for now we'll just do SU(2)

# will eventually need to replace with boilerplate code for random SU(N) generators
const σ₁ = Complex{Float64}[0 1; 1 0]
const σ₂ = Complex{Float64}[0 -im; im 0]
const σ₃ = Complex{Float64}[1 0; 0 -1]

const σ⃗ = [σ₁, σ₂, σ₃]


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


# TODO this version is probably ridiculously inefficient
"""
    randSU2()

Creates a random SU(2) matrix, uniformly distributed
on SU(2).
"""
randSU2(ϕ::AbstractFloat, n::AbstractVector) = expm(im*ϕ*n'*σ⃗)
randSU2() = randSU2(randazimuth(), randunitvec(3))
export randSU2
