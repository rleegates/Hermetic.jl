module DensePolynomials

import Combinatorics: doublefactorial
import Base: *, +, ^, scale!, size, show, string, convert
import StaticArrays: SVector

include("./DensePolynomials/mono.jl")
include("./DensePolynomials/fastbinomial.jl")
include("./DensePolynomials/smallpolynomials.jl")

immutable ProductPoly{M,T<:Number}
    k::Int
    o::Int
    c::Vector{T}
    function ProductPoly(k::Int)
        (k, o, c) = _set_ppoly(T, M, k)
        return new(k, o, c)
    end
    function ProductPoly(k::Int, c::Vector{T})
        (k, o, c) = _set_ppoly(c, M, k)
        return new(k, o, c)
    end
end

function _set_ppoly{T<:Number}(::Type{T}, m::Int, k::Int)
    na = fast_binomial(m+k, k)
    c = zeros(T,na)
    (k, na, c)
end

function _set_ppoly{T<:Number}(c::Vector{T}, m::Int, k::Int)
    na = fast_binomial(m+k, k)
    (k, na, c)
end

function string{M,T<:Number}(p::ProductPoly{M,T})
     str = "Dimension: $(M), Order: $(p.k)\nP(x) = "
     f = zeros(Int,M)
     f[1] = p.k
     i = p.o
     str *= "\n   $(round(p.c[i],2))\t* x^$f"
     while i > 1
          unsafe_mono_last_grlex!(f,M)
          i -= 1
          str *= "\n + $(round(p.c[i],2))\t* x^$f"
     end
     return str
end

function show(io::IO, p::ProductPoly)
     print(io,string(p))
end

function polynomial_sort!{F <: Real}(c::Array{F, 1}, e::Array{Int,1})
     i = sortperm(e)
     c[:] = c[i]
     e[:] = e[i]
end

function (+){M,T<:Number}(p::ProductPoly{M,T},q::ProductPoly{M,T})
    if p.k >= q.k
        k = p.k
        _k = q.k
        c = p.c
        _c = q.c
        o = p.o
        _o = q.o
    else
        k = q.k
        _k = p.k
        c = q.c
        _c = p.c
        o = q.o
        _o = p.o
    end
     c_res = zeros(T,o)
     @simd for i = 1:_o
          @inbounds c_res[i] = _c[i] + c[i]
     end
     @simd for i = _o+1:o
          @inbounds c_res[i] = c[i]
     end
     return ProductPoly{M,T}(k,c_res)
end

import Base.deepcopy
deepcopy{M,T<:Number}(p::ProductPoly{M,T}) = ProductPoly{M,T}(p.k, copy(p.c))

function (+){M,T<:Number}(p::ProductPoly{M,T},c::T)
	pp = deepcopy(p)
	pp.c[1] += c
	return pp 
end

function (+){M,T<:Number}(c::T,p::ProductPoly{M,T})
	return p+c
end


function (*){M,T<:Number}(p::ProductPoly{M,T},q::ProductPoly{M,T})
    k = p.k+q.k
    m = M
    nm = fast_binomial(m+k,k) 
    sero = zero(T)
    f  = Vector{Int}(m)
    f1 = zeros(Int,m)
    f2 = zeros(Int,m)
    c  = zeros(T,nm)
    @inbounds for j = 1:q.o
        for i = 1:p.o
            for k = 1:m
                f[k] = f1[k] + f2[k]
            end
            idx = mono_rank_grlex(m, f)
            c[idx] += p.c[i] * q.c[j]
            unsafe_mono_next_grlex!(f2, m)
        end
        unsafe_mono_next_grlex!(f1, m)
        fill!(f2,sero)
    end
    return ProductPoly{M,T}(sum(f),c)
end

function (*){M,T<:Number}(p::ProductPoly{M,T}, c::T)
     p_res = ProductPoly{M,T}(p.k)
     @simd for i = 1:p.o
          @inbounds p_res.c[i] = c*p.c[i]
     end
     return p_res
end

function (*){M,T<:Number}(c::T, p::ProductPoly{M,T})
	return p*c
end

function (^){M,T<:Number}(p::ProductPoly{M,T}, j::Integer)
    if j == 0
        return ProductPoly{M,T}(0)
    elseif j == 1
        return deepcopy(p)
    else
        return foldl(*,[p for i = 1:j])
    end
    return p
end

function dif{M,T<:Float64}(p::ProductPoly{M,T},n::Int)
	k = p.k
	nαs_old = fast_binomial(M+k,M)
	nαs_new = nαs_old - fast_binomial(M+k-1,k)
	c = zeros(T,nαs_new)
	f = zeros(Int,M)
	f[1] = k
	for i = 0:nαs_old-1
		mono_rank = nαs_old - i
		if f[n] > 0  
			f[n] -= 1
			rank = mono_rank_grlex(M,f)
			c[rank] += f[n]*p.c[mono_rank]
			f[n] += 1
		end
		unsafe_mono_last_grlex!(f,M)
	end
	#k = p.k
	#f = zeros(Int, M)
    #f[1] = k
    #o = fast_binomial(M+k,k)
    #_o = o-fast_binomial(M+k-1,k)
    #u = zeros(T,_o)
    #mono_rank = o
    #mm = M-1
    #r = k
    #tz_offs = zero(Int)
    #for r_rev = 0:(k-1)
    #    r = k - r_rev
    #    tz_offs = zero(Int)
    #    for d_rev = 0:mm
    #        d = M - d_rev
    #        m_offs = monomial_offset(r,M,n)
    #        offs = m_offs + tz_offs
    #        for _ = 1:fast_binomial(d+r-2,d-1)
    #           	mind = mono_rank-offs 
    #           	if f[n] > 0            
    #           	u[mind] += f[n]*p.c[mono_rank]
    #           end
    #           	unsafe_mono_last_grlex!(f,M)
    #           	mono_rank -= 1
    #        end
    #        if d > 1; tz_offs += trailing_zero_offset(r,M,n+1); end            
    #    end
    #end
    return ProductPoly{M,T}(k-1,c)
end

function generate_evaluate(max_m::Int,max_k::Int)
	str = "function fast_evaluate{M,T<:Number}(p::ProductPoly{M,T},x::AbstractVector)\n"
	str *= "c = p.c \n"
	str *= "k = p.k \n"
	str *= "if k == 1 \n"
	str *="	return horner_1(Val{M},c,x)\n"
	for k = 2:max_k
		str *= "elseif k == $k \n"
		str*="	return horner_$(k)(Val{M},c,x)\n"
		end
	str *= "end \n"
	str *= "end \n"
	return str
end

function generate_horner_rule(m::Int, k::Int)
    f = zeros(Int, m)
    f[1] = k
    o = fast_binomial(m+k,k)
    _o = o-fast_binomial(m+k-1,k)
    u = Vector{String}(_o)
    for i = 1:_o
    	u[i] = "c[$i]"
    end
    mono_rank = o
    mm = m-1
     r = k
     tz_offs = zero(Int)
     for d_rev = 0:mm
          d = m - d_rev
          i = d_rev + 1
          m_offs = monomial_offset(r,m,i)
          offs = m_offs + tz_offs
          for _ = 1:fast_binomial(d+r-2,d-1)
                mind = mono_rank-offs
                #@simd for j = 1:nvals
                     #@inbounds u[j,mind] += c[mono_rank] * xstat[j][i]
                #end
                umind = u[mind]            
                u[mind] = "($umind + c[$mono_rank] * x[$i])"
                unsafe_mono_last_grlex!(f,m)
                mono_rank -= 1
          end
          if d > 1; tz_offs += trailing_zero_offset(r,m,i+1); end
     end
    for r_rev = 1:(k-1)
        r = k - r_rev
        tz_offs = zero(Int)
        for d_rev = 0:mm
            d = m - d_rev
            i = d_rev + 1
            m_offs = monomial_offset(r,m,i)
            offs = m_offs + tz_offs
            for _ = 1:fast_binomial(d+r-2,d-1)
                mind = mono_rank-offs
                     #@simd for j = 1:nvals
                     #     @inbounds u[j,mind] += u[j,mono_rank] * xstat[j][i]
                     #end
                     umind = u[mind]
                     umono = u[mono_rank]
                     u[mind] = "($umind + $umono * x[$i])"
                 unsafe_mono_last_grlex!(f,m)
                 mono_rank -= 1
            end
            if d > 1; tz_offs += trailing_zero_offset(r,m,i+1); end
        end
    end
    return u[1]
end

macro generate_evaluate_schemes(max_m::Int,max_k::Int)
	fun = parse(generate_evaluate(max_m,max_k))
	return quote
	$fun
	end
end

macro generate_horner_rules(max_m::Int,max_k::Int)
	funmat = Matrix{Function}(max_m,max_k+1)
	for m = 1:max_m
		for k = 1:max_k
			eval(parse("function horner_$(k)(m::Type{Val{$m}},c::AbstractVector,x::AbstractVector); return $(generate_horner_rule(m,k)); end"))
		end
	end
end

@generated function evaluate{M,T<:Number}(p::ProductPoly{M,T}, x::AbstractVector)
    @generate_horner_rules(5,10)
    @generate_evaluate_schemes(5,10)
    return quote
    	return $(fast_evaluate)(p,x) 
	end
end

function fill_coeffs{T<:Real}(c::Vector{T}, d::Int, r::Int, nvals::Int)
     o = fast_binomial(d+r,r)-fast_binomial(d+r-1,r)
     cc = Matrix{T}(nvals, o)
     for j = 1:o
          @inbounds cj = c[j]
          @simd for i = 1:nvals
                @inbounds cc[i,j] = cj
          end
     end
     return cc
end

function fill_coeffs{M,T<:Number}(::Type{ProductPoly{M,T}}, c::Vector{T}, d::Int, r::Int, nvals::Int)
     o = fast_binomial(d+r,r)-fast_binomial(d+r-1,r)
     cc = Matrix{ProductPoly{M,T}}(nvals, o)
     for j = 1:o
          @inbounds cj = c[j]
          @simd for i = 1:nvals
                p = ProductPoly{M,T}(d,0)
                p.c[1] = cj
                @inbounds cc[i,j] = p
          end
     end
     return cc
end

function polynomial_value{T<:Number}(m::Int, o::Int,c::Array{T, 1},x::Vector{AbstractVector})
     nvals = length(x)
     f = zeros(Int,m)
     p = zeros(T, nvals)
     uno = one(T)
     v = ones(T,nvals)
     @simd for i = eachindex(p)
          @inbounds p[i] += c[1]
     end
     for j = 2:o
        unsafe_mono_next_grlex!(f,m)
        cj = c[j]
        if abs(cj) < 1e-16
            fill!(v,uno)
            for k = 1:nvals
                for i = 1:N
                     @inbounds v[k] *= x[k][i]^f[i]
                end
            end
            @simd for i = eachindex(p)
                @inbounds p[i] += cj*v[i]
            end
        end 
     end
     return p
end

function polynomial_value_horner_rule{T<:Number, I<:Int, N}(m::I, k::I, o::I, c::Array{T, 1}, xstat::Vector{SVector{N,T}})
    nvals = length(xstat)
    f = zeros(I, m)
    f[1] = k
    u = fill_coeffs(c,m,k,nvals)
    mono_rank = o
    mm = m-1
     r = k
     tz_offs = zero(I)
     for d_rev = 0:mm
          d = m - d_rev
          i = d_rev + 1
          m_offs = monomial_offset(r,m,i)
          offs = m_offs + tz_offs
          for _ = 1:fast_binomial(d+r-2,d-1)
                mind = mono_rank-offs
                @simd for j = 1:nvals
                     @inbounds u[j,mind] += c[mono_rank] * xstat[j][i]
                end
                unsafe_mono_last_grlex!(f,m)
                mono_rank -= 1
          end
          if d > 1; tz_offs += trailing_zero_offset(r,m,i+1); end
     end
    for r_rev = 1:(k-1)
        r = k - r_rev
        tz_offs = zero(I)
        for d_rev = 0:mm
            d = m - d_rev
            i = d_rev + 1
            m_offs = monomial_offset(r,m,i)
            offs = m_offs + tz_offs
            for _ = 1:fast_binomial(d+r-2,d-1)
                mind = mono_rank-offs
                     @simd for j = 1:nvals
                          @inbounds u[j,mind] += u[j,mono_rank] * xstat[j][i]
                     end
                 unsafe_mono_last_grlex!(f,m)
                 mono_rank -= 1
            end
            if d > 1; tz_offs += trailing_zero_offset(r,m,i+1); end
        end
    end
    return u[:,1]
end

function polynomial_value_horner_rule{T<:Number, I<:Int, N}(m::I, k::I, o::I,c::Array{T, 1},xstat::Vector{SVector{N,ProductPoly{T,I}}})
     nvals = length(xstat)
     f = zeros(I, m)
     f[1] = k
     u = fill_coeffs(ProductPoly{T,I},c,m,k,nvals)
     mono_rank = o
     mm = m-1
     r = k
     tz_offs = zero(I)
     for d_rev = 0:mm
          d = m - d_rev
          i = d_rev + 1
          m_offs = monomial_offset(r,m,i)
          offs = m_offs + tz_offs
          for _ = 1:fast_binomial(d+r-2,d-1)
                mind = mono_rank-offs
                @simd for j = 1:nvals
                     fac = xstat[j][i] * c[mono_rank]
                     uj = u[j,mind]
                     @inbounds u[j,mind] = uj + fac
                end
                unsafe_mono_last_grlex!(f,m)
                mono_rank -= 1
          end
          if d > 1; tz_offs += trailing_zero_offset(r,m,i+1); end
     end
     for r_rev = 1:(k-1)
          r = k - r_rev
          tz_offs = zero(I)
          for d_rev = 0:mm
                d = m - d_rev
                i = d_rev + 1
                m_offs = monomial_offset(r,m,i)
                offs = m_offs + tz_offs
                for _ = 1:fast_binomial(d+r-2,d-1)
                     mind = mono_rank-offs
                     @simd for j = 1:nvals
                          @inbounds u[j,mind] += u[j,mono_rank] * xstat[j][i]
                     end
                     unsafe_mono_last_grlex!(f,m)
                     mono_rank -= 1
                end
                if d > 1; tz_offs += trailing_zero_offset(r,m,i+1); end
          end
     end
     return u[:,1]
end

end #module DensePolynomials

import DensePolynomials: ProductPoly, evaluate
import MultiPoly: MPoly, generators
using BenchmarkTools
import StaticArrays: SVector

#p = ProductPoly{Float64,Int}(2,1)
#fill!(p.c,1.)
#x₀ = SVector{2,Float64}([1.,2.])
#@benchmark evaluate(p,x₀)

#I = Int64
#
#x, y, z = generators(MPoly{Float64}, :x, :y, :z);
#for i = 1:0
#info("polynom of order $i")
#p1 = sum([(x+y+z)^i for i = 0:i])
#foreach(x->p1.terms[x]=1.,keys(p1.terms))
#p2 = ProductPoly{Float64,I}(I(3),I(i))
#fill!(p2.c,1.)
#info("MPoly add")
#display(@benchmark $p1+$p1)
#info("DensePolynomials add")
#display(@benchmark $p2+$p2)
#info("MPoly mul")
#display(@benchmark $p1*$p1)
#info("DensePolynomials mul")
#display(@benchmark $p2*$p2)
#
#x₀_svec = Vector{SVector{3,Float64}}(1)
#x₀_svec[1] = [1.,2.,3]
#x₀ = [1.,2.,3]
#info("MPoly evaluate @ Float64")
#display(@benchmark MultiPoly.evaluate($p1,$(x₀[1]),$(x₀[2]),$(x₀[3])))
#info("DensePolynomials evaluate @ Float64")
#display(@benchmark evaluate($p2,$x₀_svec))
#end
#
#p = ProductPoly{Float64,Int}(2,2)
#p.c[1:6] = 1:6
#x₀_svec = Vector{SVector{2,Float64}}(1)
#x₀_svec[1] = [1.,2.]
#evaluate(p,x₀_svec)
#
#m = 4
#k = 4
#o = binomial(m+k,k)
#c = collect(1:o)
#x = [1.,2.]
#import DensePolynomials: fast_binomial, generate_functions, generate_horner_rule
#generate_functions(3,3)
#
#
#pp = ProductPoly{Float64,Int}(2,2)
#pp.c[1:6] = 1:6
#literals = [:a,:b,:c,:d,:e,:f,:g,:h,:i,:j]
#for m = 1:5
#	for k = 2:5
#		x = generators(MPoly{Float64}, literals[1:m]...);	
#		sumpoly = foldl(+,1,x)
#		p_mpoly = foldl(*,[sumpoly for i = 1:k])
#		foreach(x->p_mpoly.terms[x]=1.,keys(p_mpoly.terms))
#		p_dp = ProductPoly{m,Float64}(k)
#		fill!(p_dp.c,1.)
#		xx = [sumpoly for i = 1:m]
#		prpoly = ProductPoly{m,Float64}(1)
#		xxx = [prpoly for i = 1:m]
#		#nvals = ntuple(x->vals[i],m)
#		info("m = $m, k = $k")
#		info("Mpoly eval @add")
#		display(@benchmark $p_mpoly+$p_mpoly)
#		info("DensePolynomials add")
#		display(@benchmark $p_dp+$p_dp)
#		info("Mpoly mul @Float64")
#		display(@benchmark $p_mpoly*$p_mpoly)
#		info("DensePolynomials mul @Float64")
#		display(@benchmark $p_dp*$p_dp)
#		info("Mpoly eval @Float64")
#		display(@benchmark MultiPoly.evaluate($p_mpoly,$(rand(m))...))
#		info("DensePolynomials eval @Float64")
#		display(@benchmark evaluate($p_dp, $(rand(m))))
#		info("Mpoly eval @poly")
#		display(@benchmark MultiPoly.evaluate($p_mpoly,$xx...))
#		info("DensePolynomials eval @poly")
#		display(@benchmark evaluate($p_dp, $xxx))
#		println()
#	end
#end

m = 2
k = 1
p_dp = ProductPoly{m,Float64}(k)

#f@profile for i = 1:10000
#fevaluate(p,x₀)
#fend
#fusing ProfileView
#fProfileView.view()
#x₀ = SVector{2,ProductPoly{Float64,Int}}
#x₀_svec[1] = [p,p]
#pp = ProductPoly{2,Float64}(3)
#prpoly = ProductPoly{2,Float64}(1)
#xxx = [prpoly for i = 1:2]
#using ProfileView
#Profile.clear()
#@profile for i = 1:1000
#	evaluate(pp,xxx)
#end
#ProfileView.view()
#
#pp = ProductPoly{2,Float64}(3)
#prpoly = ProductPoly{2,Float64}(1)
#xxx = [rand() for i = 1:2]

#
#@benchmark evaluate(pp,x₀_svec)
#
#display(ProfileView.view())	
#Profile.clear()
#@profile for i = 1:1000
#	DensePolynomials.horner_2_2(pp.c,[p,p])
#end
#ProfileView.view()

#horner = generate_horner_rule(m,k)


#info("MPoly evaluate @ MPoly")
#display(@benchmark MultiPoly.evaluate($p1,$p1,$p1,$p1))
#info("DensePolynomials evaluate @ DensePolynomial")
#
#x₀_poly = Vector{SVector{3,ProductPoly{Float64,I}}}(1)
#x₀_poly[1] = [p2,p2,p2]
#display(@benchmark evaluate($p2,$x₀_poly))

#Profile.clear()
#@profile for i = 1:100000; evaluate(p2,x₀_svec); end
#using ProfileView
#ProfileView.view()

#Profile.clear()
#@profile for i = 1:100000; sortperm(rand(Int,5)); end
#using ProfileView
#ProfileView.view(C=true)

