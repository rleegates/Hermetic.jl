function do_plus{T<:Number}(a::T,b::T,c::T,d::T,e::T,f::T,g::T,h::T,i::T,j::T,k::T,l::T,m::T,n::T,o::T,p::T,q::T,r::T,s::T,t::T,u::T,v::T,w::T,x::T,y::T,z::T)
	return a+b+c+d+e+f+g+h+i+j+k+l+m+n+o+p+q+r+s+t+u+v+w+x+y+z
end

function do_plus2{T<:Number}(a::T,b::T,c::T,d::T,e::T,f::T,g::T,h::T,i::T,j::T,k::T,l::T,m::T,n::T,o::T,p::T,q::T,r::T,s::T,t::T,u::T,v::T,w::T,x::T,y::T,z::T)
	return (((((((((((((((((((((((((a+b)+c)+d)+e)+f)+g)+h)+i)+j)+k)+l)+m)+n)+o)+p)+q)+r)+s)+t)+u)+v)+w)+x)+y)+z)
end

function test_plus()
	for i = 1:1_000_000
		do_plus(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)
	end
	return nothing
end

function test_plus2()
	for i = 1:1_000_000
		do_plus2(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)
	end
	return nothing
end

versioninfo()
# Julia Version 0.5.1
# Commit 6445c82 (2017-03-05 13:25 UTC)
# Platform Info:
#   OS: macOS (x86_64-apple-darwin13.4.0)
#   CPU: Intel(R) Core(TM) i7-4850HQ CPU @ 2.30GHz
#   WORD_SIZE: 64
#   BLAS: libopenblas (USE64BITINT DYNAMIC_ARCH NO_AFFINITY Haswell)
#   LAPACK: libopenblas64_
#   LIBM: libopenlibm
#   LLVM: libLLVM-3.7.1 (ORCJIT, haswell)

test_plus()
test_plus2()

using BenchmarkTools

@benchmark test_plus()
# julia> @benchmark test_plus()
# BenchmarkTools.Trial:
#   memory estimate:  0 bytes
#   allocs estimate:  0
#   --------------
#   minimum time:     71.984 ms (0.00% GC)
#   median time:      73.310 ms (0.00% GC)
#   mean time:        74.043 ms (0.00% GC)
#   maximum time:     83.317 ms (0.00% GC)
#   --------------
#   samples:          68
#   evals/sample:     1


@benchmark test_plus2()
# julia> @benchmark test_plus2()
# BenchmarkTools.Trial:
#   memory estimate:  0 bytes
#   allocs estimate:  0
#   --------------
#   minimum time:     7.079 ms (0.00% GC)
#   median time:      7.252 ms (0.00% GC)
#   mean time:        7.395 ms (0.00% GC)
#   maximum time:     11.419 ms (0.00% GC)
#   --------------
#   samples:          676
#   evals/sample:     1



# Profile.clear()
# @profile test_plus()
# using ProfileView
# ProfileView.view(C=true)
