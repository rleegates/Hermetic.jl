function fast_binomial(n::Int,k::Int)
     if 0 <= k <= n
        ntok = 1
        ktok = 1
        for t in 1 : (min(k, n - k)+1)
            ntok *= n
            ktok *= t
            n -= 1
        end
        #println(ntok)
        #println(ktok)
        return div(ntok,ktok)
    else
        return 0
    end
end

function fast_binomial_2( n::Int, k::Int )
    if 0 <= k <= n
        nums = 1
        dens = 1
        for i = 1:k  
            nums *= ( n - k + i )
            dens *= i
        end
        return div(nums,dens)
    else
        return 0
    end
end

function fast_binomial_3(n::Int,k::Int)
    if 0 <= k <= n
        t = 1
        if k<n-k
            for i_rev = 0:(n-k)
                i = n-i_rev
                t = div(t*i,(n-i+1))
            end
        else
            for i_rev = 0:k
                i = n-i_rev
                t = div(t*i,(n-i+1))
            end
        end
        return t
    else
        return 0
    end
end

function recursive_binomal(n,k)
  n >= k || return 0              #short circuit base cases
  n == 1 && return 1
  k == 0 && return 1
 
  (n * recursive_binomal(n - 1, k - 1)) ÷ k  #recursive call
end

function binomial_matrix(n_max::Int, k_max::Int)
    binmat = Matrix{Int}(n_max+1,k_max+1)
    for i = 0:n_max
        for j = 0:k_max
            binmat[i+1,j+1] = binomial(i,j)
        end
    end
    return binmat
end

begin
    local binmat = binomial_matrix(50,50)
    function fastest_binomial(n::Int,k::Int)
        return binmat[n+1,k+1]
    end
end


function count_time{F<:Function}(binomi::F,n::Int,k::Int)
    then = now()
    for i = 1:1000000
        binomi(n,k)
    end
    return convert(Int,(now() - then))
end

count_time_binomial(n::Int,k::Int) = count_time(binomial,n,k)
count_time_fast_binomial(n::Int,k::Int) = count_time(fast_binomial,n,k)
count_time_fast_binomial_2(n::Int,k::Int) = count_time(fast_binomial_2,n,k)
count_time_fast_binomial_2(n::Int,k::Int) = count_time(fast_binomial_2,n,k)
count_time_fast_binomial_3(n::Int,k::Int) = count_time(fast_binomial_3,n,k)
count_time_fastest_binomial(n::Int,k::Int) = count_time(fastest_binomial,n,k)
count_time_recursive_binomial(n::Int,k::Int) = count_time(recursive_binomal,n,k)

function generate_pascal{F<:Function}(binomi::F,n::Int,k::Int)
    retmat = Matrix{Int}(n+1,k+1)
    for i = 0:n
        for j = 0:k
            retmat[i+1,j+1] = binomi(i,j)
        end
    end
    return retmat
end

function generate_n_minus_k(n::Int,k::Int)
    retmat = Matrix{Int}(n+1,k+1)
    for i = 0:n
        for j = 0:k
            retmat[i+1,j+1] = i-j
        end
    end
    return retmat
end

function generate_k(n::Int,k::Int)
    retmat = Matrix{Int}(n+1,k+1)
    for i = 0:n
        for j = 0:k
            retmat[i+1,j+1] = j
        end
    end
    return retmat
end

N = 10
using BenchmarkTools
display(@benchmark generate_pascal(binomial,N,N))
display(@benchmark generate_pascal(fast_binomial,N,N))
display(@benchmark generate_pascal(fast_binomial_2,N,N))
display(@benchmark generate_pascal(fast_binomial_3,N,N))
display(@benchmark generate_pascal(fastest_binomial,N,N))
display(@benchmark generate_pascal(recursive_binomal,N,N))
display(generate_pascal(count_time_binomial,N,N))
display(generate_pascal(count_time_fast_binomial,N,N))
display(generate_pascal(count_time_fast_binomial_2,N,N))
display(generate_pascal(count_time_fastest_binomial,N,N))
display(generate_pascal(count_time_fast_binomial_3,N,N))
display(generate_pascal(count_time_recursive_binomial,N,N))