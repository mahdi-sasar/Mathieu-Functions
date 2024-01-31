using LinearAlgebra

#Characteristic exponent calculator for even indices:
function CharA(n::Int64, q::Float64)
    if n%2 == 0
        m = floor(Int, (n/2) + 1)
        Nmax = 2*m
        B = zeros(Nmax, Nmax)
        B[1, 1] = 0.
        B[1, 2] = (-1.)*sqrt(2.)*q
        B[2, 1] = (-1.)*sqrt(2.)*q
        B[2, 2] = 4.
        B[2, 3] = (-1.)*q
        B[3, 2] = (-1.)*q
        for i in 3:Nmax
            for j in 3:Nmax
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
        Nmax = 2*m
        B = zeros(Nmax, Nmax)
        B[1, 1] = (1.)*(q + 1.)
        B[1, 2] = (-1.)*q
        B[2, 1] = (-1.)*q
        B[2, 2] = 9.
        B[2, 3] = (-1.)*q
        B[3, 2] = (-1.)*q
        for i in 3:Nmax
            for j in 3:Nmax
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
        Nmax = 2*m
        B = zeros(Nmax, Nmax)
        B[1, 1] = 4.
        B[1, 2] = -q
        B[2, 1] = (-1.)*q
        B[2, 2] = 16.
        B[2, 3] = (-1.)*q
        B[3, 2] = (-1.)*q
        for i in 3:Nmax
            for j in 3:Nmax
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
        Nmax = 2*m
        B = zeros(Nmax, Nmax)
        B[1, 1] = 1. - q
        B[1, 2] = -q
        B[2, 1] = -q
        B[2, 2] = 9.
        B[2, 3] = -q
        B[3, 2] = -q
        for i in 3:Nmax
            for j in 3:Nmax
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
