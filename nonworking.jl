using LinearAlgebra,Random, Plots, NonlinearSolve
using StaticArrays, ForwardDiff, Distributed, StaticArrays, SciMLNLSolve, SimpleNonlinearSolve, LineSearches, NLsolve


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

f = NonlinearFunction(f!)

init = ones(152);
A = ones(152);
A[6] = 0.8

f = NonlinearFunction(f!, jac = j!)

p = A

ProbN = NonlinearProblem(f,init,p)
sol = solve(ProbN,NewtonRaphson(), reltol = 1e-8,abstol = 1e-8)

