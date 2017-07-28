function polynomial_value_m_k0{T<:Number, N}(coeffs::Array{T,1}, x::Vector{SVector{N,T}})
    len = length(x)
    c = coeffs[1]
    retval = Vector{T}(len)
    @simd for i = 1:len
        @inbounds retval[i] = c 
    end
    return retval
end

function polynomial_value_m1_k1{T<:Number, N}(coeffs::Array{T,1}, x::Vector{SVector{N,T}})
    len = length(x)
    a = coeffs[1]
    b = coeffs[2]
    retval = Vector{T}(len)
    @simd for i = 1:len
        @inbounds retval[i] = a + b*x[i][1] 
    end
    return retval
end

function polynomial_value_m1_k2{T<:Number, N}(coeffs::Array{T,1}, x::Vector{SVector{N,T}})
    len = length(x)
    a = coeffs[1]
    b = coeffs[2]
    c = coeffs[3]
    retval = Vector{T}(len)
    @simd for i = 1:len
        @inbounds retval[i] = ((a + b*x[i][1]) + c*x[i][1]^2)
    end
    return retval
end

function polynomial_value_m1_k3{T<:Number, N}(coeffs::Array{T,1}, x::Vector{SVector{N,T}})
    len = length(x)
    a = coeffs[1]
    b = coeffs[2]
    c = coeffs[3]
    d = coeffs[4]
    retval = Vector{T}(len)
    @simd for i = 1:len
        @inbounds retval[i] = (((a + b*x[i][1]) + c*x[i][1]^2) + d*x[i][1]^3)
    end
    return retval
end

function polynomial_value_m2_k1{T<:Number, N}(coeffs::Array{T,1}, x::Vector{SVector{N,T}})
    len = length(x)
    a = coeffs[1]
    b = coeffs[2]
    c = coeffs[3]
    retval = Vector{T}(len)
    @simd for i = 1:len
        @inbounds retval[i] = ((a + b*x[i][2]) + c*x[i][1])
    end
    return retval
end

function polynomial_value_m2_k2{T<:Number, N}(coeffs::Array{T,1}, x::Vector{SVector{N,T}})
    len = length(x)
    a = coeffs[1]
    b = coeffs[2]
    c = coeffs[3]
    d = coeffs[4]
    e = coeffs[5]
    f = coeffs[6]
    retval = Vector{T}(len)
    @simd for i = 1:len
        @inbounds retval[i] = (((((a + b*x[i][2]) + c*x[i][1]) + d*x[i][2]^2) + e*x[i][1]*x[i][2]) + f*x[i][1]^2)
    end
    return retval
end

function polynomial_value_m2_k3{T<:Number, N}(coeffs::Array{T,1}, x::Vector{SVector{N,T}})
    len = length(x)
    a = coeffs[1]
    b = coeffs[2]
    c = coeffs[3]
    d = coeffs[4]
    e = coeffs[5]
    f = coeffs[6]
    g = coeffs[7]
    h = coeffs[8]
    i = coeffs[9]
    j = coeffs[10]
    retval = Vector{T}(len)
    @simd for i = 1:len
        @inbounds retval[i] = (((((((((a + b*x[k][2]) + c*x[k][1]) + d*x[k][2]^2) + e*x[k][1]*x[k][2]) + f*x[k][1]^2) + g*x[k][2]^3) + h*x[k][1]*x[k][2]^2) + i*x[k][1]^2*x[k][2]) + j*x[k][1]^3)
    end
    return retval
end

function polynomial_value_m3_k1{T<:Number, N}(coeffs::Array{T,1}, x::Vector{SVector{N,T}})
    len = length(x)
    a = coeffs[1]
    b = coeffs[2]
    c = coeffs[3]
    d = coeffs[4]
    retval = Vector{T}(len)
    @simd for i = 1:len
        @inbounds retval[i] = a + b*x[i][3] + c*x[i][2] + d*x[i][1] 
    end
    return retval
end

function polynomial_value_m3_k2{T<:Number, N}(coeffs::Array{T,1}, x::Vector{SVector{N,T}})
    len = length(x)
    a = coeffs[1]
    b = coeffs[2]
    c = coeffs[3]
    d = coeffs[4]
    e = coeffs[5]
    f = coeffs[6]
    g = coeffs[7]
    h = coeffs[8]
    i = coeffs[9]
    j = coeffs[10]
    retval = Vector{T}(len)
    @simd for k = 1:len
        @inbounds retval[k] = (((((((((a + b*x[k][3]) + c*x[k][2]) + d*x[k][1]) + e*x[k][3]^2) + f*x[k][2]*x[k][3]) + g*x[k][2]^2) + h*x[k][1]*x[k][3]) + i*x[k][1]*x[k][2]) + j*x[k][1]^2)
    end
    return retval
end

function polynomial_value_m3_k3{T<:Number, N}(coeffs::Array{T,1}, x::Vector{SVector{N,T}})
    len = length(x)
    @inbounds begin
        c1 = coeffs[1]
        c2 = coeffs[2]
        c3 = coeffs[3]
        c4 = coeffs[4]
        c5 = coeffs[5]
        c6 = coeffs[6]
        c7 = coeffs[7]
        c8 = coeffs[8]
        c9 = coeffs[9]
        c10 = coeffs[10]
        c11 = coeffs[11]
        c12 = coeffs[12]
        c13 = coeffs[13]
        c14 = coeffs[14]
        c15 = coeffs[15]
        c16 = coeffs[16]
        c17 = coeffs[17]
        c18 = coeffs[18]
        c19 = coeffs[19]
        c20 = coeffs[20]    
        retval = Vector{T}(len)
        for k = 1:len
            retval[k] = (((((((((((((((((((c1 + c2*x[k][3]) + c3*x[k][2]) + c4*x[k][1]) + c5*x[k][3]^2) + c6*x[k][2]*x[k][3]) + c7*x[k][2]^2) + c8*x[k][1]*x[k][3]) + c9*x[k][1]*x[k][2]) + c10*x[k][1]^2) + c11*x[k][3]^3) + c12*x[k][2]*x[k][3]^2) + c13*x[k][2]^2*x[k][3]) + c14*x[k][2]^3) + c15*x[k][1]*x[k][3]^2) + c16*x[k][1]*x[k][2]*x[k][3]) + c17*x[k][1]*x[k][2]^2) + c18*x[k][1]^2*x[k][3]) + c19*x[k][1]^2*x[k][2]) + c20*x[k][1]^3)
        end
    end
    return retval
end