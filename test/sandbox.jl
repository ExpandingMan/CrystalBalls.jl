using CrystalBalls
using StaticArrays
using BenchmarkTools

lat = LatticeToroidal(2, 1.0, 3)

x = SVector(0, 0)

# U = GaugeFieldSU2(lat)
U = rand(GaugeFieldSU2, lat)
