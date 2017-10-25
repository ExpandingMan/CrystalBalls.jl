using CrystalBalls
using StaticArrays

lat = LatticeToroidal(1, 1.0, 3)

U = GaugeFieldSU2(lat)
# U = rand(GaugeFieldSU2, lat)
