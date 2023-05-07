using LinearAlgebra, NonlinearSolve
using ForwardDiff, SciMLNLSolve, LineSearches


function problem(x,A)
    return x .^ 2 - A
end

function problemJacobian(x,A)
    return diagm(2 .* x)
end

function f!(F,u,p)
    F[1:152] = problem(u,p)
end

function j!(J,u,p)
    J[1:152,1:152] = problemJacobian(u,p)
end

f = NonlinearFunction(f!,jac = j!)

init = ones(Complex{Float64},152);
A = ones(152);
A[6] = 0.8

f = NonlinearFunction(f!, jac = j!)

p = A

ProbN = NonlinearProblem(f,init,p)
sol = solve(ProbN,NLSolveJL(linesearch = HagerZhang(),method = :newton), reltol = 1e-8,abstol = 1e-8)


print(sol.u)
eltype(init)<: Number

typeof(init)<: Vector{T} where T <: Number

