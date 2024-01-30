using LinearAlgebra

#Characteristic exponent calculator for even indices:
function CharA(n::Int64, q::Float64)
    if n%2 == 0
        m = floor(Int, (n/2) + 1)
        # Nmax is an arbitrary definition at this point.
        # Nmax = max(50+m, 10*sqrt(abs(q)))  
        Nmax = 2*m
        NmaxInt = floor(Int, Nmax)
        B = zeros(NmaxInt, NmaxInt)
        B[1, 1] = 0.
        B[1, 2] = (-1.)*sqrt(2.)*q
        B[2, 1] = (-1.)*sqrt(2.)*q
        B[2, 2] = 4.
        B[2, 3] = (-1.)*q
        B[3, 2] = (-1.)*q
        for i in 3:NmaxInt
            for j in 3:NmaxInt
                if i == j
                    B[i, j] = (4.)*(i-1)^2
                elseif (i-j == 1 || j-i == 1)
                    B[i, j] = (-1.)*q
                end
            end
        end
        c = eigvals(B)
        return c[m]
    elseif n%2 == 1
        m = floor(Int, (n+1)/2.)
        # Nmax = max(50+m, 10*sqrt(abs(q)))  
        Nmax = 2*m
        NmaxInt = floor(Int, Nmax)
        B = zeros(NmaxInt, NmaxInt)
        B[1, 1] = (1.)*(q + 1.)
        B[1, 2] = (-1.)*q
        B[2, 1] = (-1.)*q
        B[2, 2] = 9.
        B[2, 3] = (-1.)*q
        B[3, 2] = (-1.)*q
        for i in 3:NmaxInt
            for j in 3:NmaxInt
                if i == j
                    B[i, j] = (2*i - 1.)^2
                elseif (i-j == 1 || j-i == 1)
                    B[i, j] = (-1.)*q
                end
            end
        end
        c = eigvals(B)
        return c[m]
    end
end

#Characteristic exponent calculator for odd indices:
function CharB(n::Int64, q::Float64)
    if (n == 0)
        return 0.
    end
    if (n%2 == 0 && n != 0)
        m = floor(Int, (n/2))
        # Nmax = max(50+m, 10*sqrt(abs(q)))   
        Nmax = 2*m
        NmaxInt = floor(Int, Nmax)
        B = zeros(NmaxInt, NmaxInt)
        B[1, 1] = 4.
        B[1, 2] = -q
        B[2, 1] = (-1.)*q
        B[2, 2] = 16.
        B[2, 3] = (-1.)*q
        B[3, 2] = (-1.)*q
        for i in 3:NmaxInt
            for j in 3:NmaxInt
                if i == j
                    B[i, j] = (4.)*(i*i)
                elseif (i-j == 1 || j-i == 1)
                    B[i, j] = (-1.)*q
                end
            end
        end
        c = eigvals(B)
        return c[m]
    elseif n%2 == 1
        m = floor(Int, (n+1)/2.)
        # Nmax = max(50+m, 10*sqrt(abs(q)))
        Nmax = 2*m
        NmaxInt = floor(Int, Nmax)
        B = zeros(NmaxInt, NmaxInt)
        B[1, 1] = 1. - q
        B[1, 2] = -q
        B[2, 1] = -q
        B[2, 2] = 9.
        B[2, 3] = -q
        B[3, 2] = -q
        for i in 3:NmaxInt
            for j in 3:NmaxInt
                if i == j
                    B[i, j] = (2*i - 1.)*(2*i - 1.)
                elseif (i-j == 1 || j-i == 1)
                    B[i, j] = -q
                end
            end
        end
        c = eigvals(B)
        return c[m]
    end
end


#Mathieu cosine function:
function ce(m::Int64, q::Float64, x::AbstractVector{Float64})
    accuracy = 0.00001
    if m%2 == 0
        a = CharA(m, q)
        A = [0., 1., 0., a/q, 0., ((a-4.)*a/(q*q))-2., 0.]
        i = 8
        #Normalization as each new coefficient is computed:
        N = 2*A[2]^2 + sum(A[4:end].^2)
        A = A ./N
        while (abs(A[i-2]/A[2])>accuracy)
            push!(A, ((a-(i-4.)*(i-4.))/q)*A[i-2]- A[i-4])
            push!(A, 0.)
            N = 2*A[2]^2 + sum(A[4:end].^2)
            A = A ./N
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
        A = A ./N
        while (abs(A[i-2]/A[1])>accuracy)
            push!(A, ((a-(i-2)*(i-2))/q)*A[i-2] - A[i-4])
            push!(A, 0.)
            N = sum(A[1:end].^2)
            A = A ./N
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
        B = B ./N
        while (abs(B[i-2]/B[2])>accuracy)
            push!(B, ((a-(i-2.)*(i-2.))/q)*B[i-2]- B[i-4])
            push!(B, 0.)
            N = sum(B[2:end].^2)
            B = B ./N
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
        B = B ./N
        while (abs(B[i-2]/B[1])>accuracy)
            push!(B, ((a-(i-2)*(i-2))/q)*B[i-2] - B[i-4])
            push!(B, 0.)
            N = sum(B[1:end].^2)
            B = B ./N
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
