include("../src/Hermetic.jl")


import Hermetic: ProductPoly, setcoef!, polyval, unsafe_mono_next_grlex!, unsafe_mono_next_grevlex!, mono_unrank_grlex!, mono_rank_grlex, unsafe_mono_last_grlex!, polynomial_value_horner_rule
import StaticArrays: SVector
using MultiPoly
using BenchmarkTools

x1, y1, z1 = generators(MPoly{Float64}, :x, :y, :z);
p1 = sum([(x1+y1+z1)^i for i = 0:5])
foreach(x->p1.terms[x]=1.,keys(p1.terms))
x2,y2,z2,c = Hermetic.vars(ProductPoly,3)
p2 = sum([(x2+y2+z2)^i for i = 0:5])
fill!(p2.c,1.)

info("MPoly multiplicate")
display(@benchmark (p1*p1))
info("Hermetic multiplicate")
display(@benchmark (p2*p2))
info("MPoly add")
display(@benchmark (p1+p1))
info("Hermetic add")
display(@benchmark (p2+p2))


x₀_svec = Vector{SVector{3,ProductPoly{Hermetic.Standard,Int,Vector{Float64},Vector{Int}}}}(1)
x₀_svec[1] = [z2^2,y2^2,x2^3]

#p = p2
#x = x₀_svec
#@code_warntype Hermetic.polynomial_value_horner_rule(p.m, p.k, p.o, p.c, p.e, x)

#polyval(polyval(p2,x₀_svec),[1. 1.])
#evaluate(evaluate(p1,y1^2,x1^3),1.,1.)
info("MPoly evaluate Poly @ Poly")
display(@benchmark evaluate(p1,z1^2,y1^2,x1^3))
info("MPoly evaluate Poly @ Poly")
display(@benchmark polyval(p2,x₀_svec))