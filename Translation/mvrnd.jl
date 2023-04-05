
function mvnrnd(μ,Σ)

    U = cholesky(Σ)
    return Random.randn(size(μ))*U + mu
end
