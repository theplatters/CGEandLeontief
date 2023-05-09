using NLsolve, LinearAlgebra, LineSearches

function f!(F,f,x)
    F[1:152] = x .^2 .- 0.8
end

function j!(J,f,x)
    J[1:152,1:152] = diagm(2 .* x)
end

x0 = ones(ComplexF64,76)

solve(f!,j!)