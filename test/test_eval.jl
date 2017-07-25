include("../src/Hermetic.jl")


import Hermetic: ProductPoly, setcoef!, polyval, unsafe_mono_next_grlex!, unsafe_mono_next_grevlex!, mono_unrank_grlex!, mono_rank_grlex, unsafe_mono_last_grlex!
import StaticArrays: SVector
using MultiPoly
using BenchmarkTools


x, y, z, t, w = generators(MPoly{Float64}, :x, :y, :z, :t, :w);
p1 = sum([(x+y+z+t+w)^i for i = 0:5])
foreach(x->p1.terms[x]=1.,keys(p1.terms))
p2 = ProductPoly(5, 5)
#setcoef!(p2,[2.,0.,0.,0.,2.,0.,1.,0.,0.,1.])
fill!(p2.c,1.)
x₀ = [1.,2.,3.,4.,5.]
x₀_svec = Vector{SVector{5,Float64}}(1)
x₀_svec[1] = [1.,2.,3.,4.,5.]
#x₀_svec[2] = [2.,3.,4.,5.,6.]
x₀_mat = [1. 2. 3. 4. 5.]

evaluate(p1, x₀...)
polyval(p2, x₀_mat)
polyval(p2, x₀_svec)


info("MPoly evaluate")
display(@benchmark evaluate(p1, $(x₀[1]), $(x₀[2]), $(x₀[3]), $(x₀[4]), $(x₀[5])))

info("Hermetic evaluate")
display(@benchmark polyval($p2, $x₀_mat))

info("Hermetic evaluate SVector")
display(@benchmark polyval($p2, $x₀_svec))


info("MPoly evaluate several points")
xxx = rand(8, 5)
xxxx = [SVector(convert(Tuple,rand(5))) for i = 1:8]
display(@benchmark begin
    for i = 1:size(xxx,1)
        evaluate($p1, xxx[i,:]...)
    end
end)

info("Hermetic evaluate SVector several points")
display(@benchmark polyval($p2, $xxxx))


x, y, z = generators(MPoly{Float64}, :x, :y, :z);
p1 = sum([(x+y+z)^i for i = 0:5])
foreach(x->p1.terms[x]=1.,keys(p1.terms))
p2 = ProductPoly(3, 5)
#setcoef!(p2,[2.,0.,0.,0.,2.,0.,1.,0.,0.,1.])
fill!(p2.c,1.)
x₀ = [1.,2.,3.]
x₀_svec = Vector{SVector{3,Float64}}(1)
x₀_svec[1] = [1.,2.,3.]
#x₀_svec[2] = [2.,3.,4.,5.,6.]
x₀_mat = [1. 2. 3.]

evaluate(p1, x₀...)
polyval(p2, x₀_mat)
polyval(p2, x₀_svec)


info("MPoly evaluate")
display(@benchmark evaluate(p1, $(x₀[1]), $(x₀[2]), $(x₀[3])))

info("Hermetic evaluate")
display(@benchmark polyval($p2, $x₀_mat))

info("Hermetic evaluate SVector")
display(@benchmark polyval($p2, $x₀_svec))


info("MPoly evaluate several points")
xxx = rand(8, 5)
xxxx = [SVector(convert(Tuple,rand(5))) for i = 1:8]
display(@benchmark begin
    for i = 1:size(xxx,1)
        evaluate($p1, xxx[i,:]...)
    end
end)

info("Hermetic evaluate SVector several points")
display(@benchmark polyval($p2, $xxxx))


