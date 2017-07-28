
function mono_rank_grlex{T <: Int}(m::T, x::Array{T, 1})

     # @assert m > 0 "The dimension `M` must be > 1"
     # @assert m == length(x) "The monimial size is incompatible with the
     # dimension of the polynomial"

     # for i = 1 : m
     #     if (x[i] < 0)
     #         throw()
     #     end
     # end

     if m==1
          return x[1]+1
     end

     nm = sum(x)

     # Convert to KSUBSET format.

     ns = nm + m - 1;
     ks = m - 1;

     xs = cumsum(x[1:m-1]+1)

     ##  Compute the rank.

     rank = 1;

     @inbounds for i = 1:ks
          tim1 = i == 1 ? 0 : xs[i-1];

          if (tim1 + 1 <= xs[i] - 1)
                for j = tim1 + 1:xs[i] - 1
                     rank += fast_binomial(ns - j, ks - i)
                end
          end
     end

     @inbounds for n = 0:nm - 1
          rank += fast_binomial(n + m - 1, n)
     end

     return rank
end

function mono_grlex!(X::Array{Int, 2}, m::Int)
     n, s = size(X)
     @assert s==m

     X[1,:] = 0

     @inbounds for j = 2:n
          x = view(X, j, :)
          X[j,:] = mono_next_grlex!(vec(X[j-1, :]), m)
     end
     return X
end

function get_inter_idx{T <: Int}(X::Array{T, 2}, ki::T)
     n, m = size(X)
     rank = sum(X, 2)
     k = maximum(rank)
      nz = Vector{Int}(n)
      for i = 1:n
          nz[i] = sum(X[i,:] .== 0)
      end
     #nz = sum(map(x -> x == 0 ? 1 : 0, X), 2) ## Number of zero in composition

     ## Interactions are those allocation with more than 1 zero. Or, non
     ## interaction are thos allocations with exactly M-1 zeros plus the case with
     ## M zeros
     idx = BitArray(n)
     for i = 1:n
          idx[i] = nz[i] >= m-1
     end

     if ki == k
          for i = 1:n
                idx[i] = true
          end
          #idx = BitArray([1 for i = 1:n])
     elseif ki > 1 & ki < k
          ## Interactions of order ki corresponds to those allocations λᵢ with
          ## more than 1 zero and sum(λᵢ) == ki
             for i = 1:n
                idx[i] = !idx[i] & rank[i] == ki | idx[i]
             end
          #idx = BitArray((!idx & (rank .== ki)) | idx)
     end
     return idx
end

function mono_unrank_grlex{T <: Int}(m::T, rank::T)
     if (m == 1)
          return [rank-1]
     end

     rank1 = 1;
     nm = -1;
     while  true
          nm = nm + 1
          r = fast_binomial(nm + m - 1, nm)
          if (rank < rank1 + r)
                break
          end
          rank1 = rank1 + r
     end

     rank2 = rank - rank1

     ks = m - 1
     ns = nm + m - 1
     nksub = fast_binomial(ns, ks)
     xs = zeros(T, ks, 1);
     j = 1;

     @inbounds for i = 1:ks
          r = fast_binomial(ns - j, ks - i)
          while (r <= rank2 && 0 < r)
                rank2 = rank2 - r
                j = j + 1
                r = fast_binomial(ns - j, ks - i)
          end
          xs[i] = j
          j +=  1
     end

     x = zeros(T, m)
     x[1] = xs[1] - 1
     @inbounds for i = 2:m - 1
          x[i] = xs[i] - xs[i-1] - 1
     end
     x[m] = ns - xs[ks];

     return x
end

function mono_unrank_grlex!{T <: Int}(x::Array{T, 1}, m::T, rank::T)
     if (m == 1)
          x[1] = rank-1
          return x
     end

     rank1 = 1;
     nm = -1;
     while  true
          nm = nm + 1
          r = fast_binomial(nm + m - 1, nm)
          if (rank < rank1 + r)
                break
          end
          rank1 = rank1 + r
     end

     rank2 = rank - rank1

     ks = m - 1
     ns = nm + m - 1
     nksub = fast_binomial(ns, ks)
     xs = zeros(T, ks, 1);
     j = 1;

     @inbounds for i = 1:ks
          r = fast_binomial(ns - j, ks - i)
          while (r <= rank2 && 0 < r)
                rank2 = rank2 - r
                j = j + 1
                r = fast_binomial(ns - j, ks - i)
          end
          xs[i] = j
          j +=  1
     end

     x[1] = xs[1] - 1
     @inbounds for i = 2:m - 1
          x[i] = xs[i] - xs[i-1] - 1
     end
     x[m] = ns - xs[ks];

     return x
end

function mono_next_grlex!{T <: Int}(x::Array{T, 1}, m::T)
     @assert m >= 0
     @assert all(x.>=0)

     i = 0
     @inbounds for j = m:-1:1
          if 0 < x[j]
                i = j
                break
          end
     end

     if i == 0
          x[m] = 1
          return x
     elseif i == 1
          t = x[1] + 1
          im1 = m
     elseif 1 < i
          t = x[i]
          im1 = i - 1
     end

     @inbounds x[i] = 0
     @inbounds x[im1] = x[im1] + 1
     @inbounds x[m] = x[m] + t - 1

     return x
end

function unsafe_mono_next_grlex!{T <: Int}(x::Array{T, 1}, m::T)
     i = 0
     @inbounds for j = m:-1:1
          if 0 < x[j]
                i = j
                break
          end
     end

     if i == 0
          @inbounds x[m] = 1
          return x
     elseif i == 1
          @inbounds t = x[1] + 1
          im1 = m
     elseif 1 < i
          @inbounds t = x[i]
          im1 = i - 1
     end

     @inbounds x[i] = 0
     @inbounds x[im1] = x[im1] + 1
     @inbounds x[m] = x[m] + t - 1

     return x
end

findfirst_reverse_uneq_zero{T,N}(collection::AbstractArray{T,N}, start_idx::Int) = findfirstop_reverse(!=,zero(T),collection,start_idx)

function findfirstop_reverse{F<:Function,T,N}(op::F, what::T, collection::AbstractArray{T,N}, start_idx::Int)
    @inbounds for j = 0:(start_idx-1)
        i = start_idx - j
     #@inbounds for i = start_idx:-1:1
        if op(what,collection[i])
            return i
        end
     end
     return 0
end

function unsafe_mono_last_grlex!{T <: Int}(x::Array{T, 1}, m::T)
     i = findfirst_reverse_uneq_zero(x,m)
     if i < m
          @inbounds x[i] -= 1
          @inbounds x[i+1] += 1
     elseif i == m
          j = findfirst_reverse_uneq_zero(x,m-1)
          @inbounds t = x[m]
          @inbounds x[m] = 0
          if j > 0
                @inbounds x[j] -= 1
                @inbounds x[j+1] = t + 1
          else
                @inbounds x[1] = t - 1
          end
     end

     return x
end


function unsafe_mono_next_grevlex!{T <: Int}(x::Array{T, 1}, m::T)
     j = 1;
     @inbounds for i = 2 : m
          if 0 < x[i]
                j = i
                break
          end
     end

     if j == 1
          @inbounds t = x[1]
          @inbounds x[1] = 0
          @inbounds x[m] = t + 1
     elseif j < m
          @inbounds x[j] = x[j] - 1
          @inbounds t = x[1] + 1
          @inbounds x[1] = 0
          @inbounds x[j-1] = x[j-1] + t;
     elseif j == m
          @inbounds t = x[1]
          @inbounds x[1] = 0
          @inbounds x[j-1] = t + 1
          @inbounds x[j] = x[j] - 1;
     end

     return x
end

function monomial_offset(r::Integer,m::Integer,i::Integer)
   offset = fast_binomial(r+m-i,r)
   return offset
end

function trailing_zero_offset(r::Integer,m::Integer,i::Integer)
  offset = fast_binomial(r+m-1-i,m-i+1)
  return offset
end