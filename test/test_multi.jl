include("../src/Hermetic.jl")


import Hermetic: ProductPoly, setcoef!, polyval, unsafe_mono_next_grlex!, unsafe_mono_next_grevlex!, mono_unrank_grlex!, mono_rank_grlex, unsafe_mono_last_grlex!, polynomial_value_horner_rule
import StaticArrays: SVector
using MultiPoly
using BenchmarkTools

x1, y1, z1 = generators(MPoly{Float64}, :x, :y, :z);
p1 = sum([(x1+y1+z1)^i for i = 0:2])
foreach(x->p1.terms[x]=1.,keys(p1.terms))
x2,y2,z2,c = Hermetic.vars(ProductPoly,3)
p2 = sum([(x2+y2+z2)^i for i = 0:2])
fill!(p2.c,1.)

info("MPoly multiplicate")
display(@benchmark (p1*p1))
info("Hermetic multiplicate")
display(@benchmark (p2*p2))
info("MPoly add")
display(@benchmark (p1+p1))
info("Hermetic add")
display(@benchmark (p2+p2))


x₀_poly = Vector{SVector{3,ProductPoly{Hermetic.Standard,Int,Vector{Float64},Vector{Int}}}}(1)
x₀_poly[1] = [p2,p2,p2]

x₀_svec = Vector{SVector{3,Float64}}(1)
x₀_svec[1] = [1.,2.,3.]

#Profile.clear()
#@profile for i = 1:100000000; polyval(p2,x₀_svec); end
#using ProfileView
#Profile.print()

#p = p2
#x = x₀_svec
#@code_warntype Hermetic.polynomial_value_horner_rule(p.m, p.k, p.o, p.c, p.e, x)

info("MPoly evaluate Poly @ Poly")
display(@benchmark evaluate($p1,$p1,$p1,$p1))
info("MPoly evaluate Poly @ Poly")
display(@benchmark polyval($p2,$x₀_poly))
