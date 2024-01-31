include("CharExponents.jl")

#Mathieu cosine function:
function ce(m::Int64, q::Float64, x::AbstractVector{Float64})
    accuracy = 0.00001
    if m%2 == 0
        a = CharA(m, q)
        A = [0., 1., 0., a/q, 0., ((a-4.)*a/(q*q))-2., 0.]
        i = 8
        #Normalization as each new coefficient is computed:
        N = 2*A[2]^2 + sum(A[4:end].^2)
        A = A ./sqrt(N)
        while (abs(A[i-2])>accuracy)
            push!(A, ((a-(i-4.)*(i-4.))/q)*A[i-2]- A[i-4])
            push!(A, 0.)
            N = 2*A[2]^2 + sum(A[4:end].^2)
            A = A ./sqrt(N)
            i += 2
        end
        t = zeros(length(x))
        for j in 2:length(A)
            t += A[j]*cos.((j-2)*x)
        end
        return t
    end
    if m%2 == 1
        a = CharA(m, q)
        A = [1., 0., (a-q-1.)/q, 0., ((a-9.)/q)*((a-q-1.)/q)-1., 0., ((a-25.)/q)*(((a-9.)/q)*((a-q-1.)/q)-1.)-(a-q-1.)/q, 0.]
        i = 9
        #Normalization as each new coefficient is computed:
        N = sum(A[1:end].^2)
        A = A ./sqrt(N)
        while (abs(A[i-2])>accuracy)
            push!(A, ((a-(i-2)*(i-2))/q)*A[i-2] - A[i-4])
            push!(A, 0.)
            N = sum(A[1:end].^2)
            A = A ./sqrt(N)
            i += 2
        end
        t = zeros(length(x))
        for j in 1:length(A)
            t += A[j]*cos.((j)*x)
        end
        return t
    end
end

#Mathieu sine function:
function se(m::Int64, q::Float64, x::AbstractVector{Float64})
    accuracy = 0.00001
    if m%2 == 0
        a = CharB(m, q)
        B = [0., 1., 0., (a-4.)/q, 0.]
        i = 6
        #Normalization as each new coefficient is computed:
        N = sum(B[2:end].^2)
        B = B ./sqrt(N)
        while (abs(B[i-2])>accuracy)
            push!(B, ((a-(i-2.)*(i-2.))/q)*B[i-2]- B[i-4])
            push!(B, 0.)
            N = sum(B[2:end].^2)
            B = B ./sqrt(N)
            i += 2
        end
        t = zeros(length(x))
        for j in 2:length(B)
            t += B[j]*sin.((j)*x)
        end
        return t
    end
    if m%2 == 1
        a = CharB(m, q)
        B = [1., 0., (a+q-1.)/q, 0.]
        i = 5
        #Normalization as each new coefficient is computed:
        N = sum(B[1:end].^2)
        B = B ./sqrt(N)
        while (abs(B[i-2])>accuracy)
            push!(B, ((a-(i-2)*(i-2))/q)*B[i-2] - B[i-4])
            push!(B, 0.)
            N = sum(B[1:end].^2)
            B = B ./sqrt(N)
            i += 2
        end
        println(B)
        t = zeros(length(x))
        for j in 1:length(B)
            t += B[j]*sin.((j)*x)
        end
        return t
    end
end
