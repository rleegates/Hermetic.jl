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
     local const binmat = binomial_matrix(50,50)
     function fast_binomial{I<:Integer}(n::I,k::I)
          return binmat[n+1,k+1]
     end
end