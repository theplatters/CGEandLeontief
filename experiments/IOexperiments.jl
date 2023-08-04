using Graphs, GraphMakie, GLMakie, MAT, LinearAlgebra, SimpleWeightedGraphs
using GraphMakie.NetworkLayout
include("../src/loaders.jl")

function getVariables(year)
    IO = data[data[:, 1].==year, :]
    temp = [8, 60, 62, 80:88...]
    IO = IO[setdiff(1:end, temp), setdiff(1:end, [1, 2, 3, 4, 5, 94])]
    IO = IO[:, setdiff(1:end, temp)]
    Ω = IO ./ sum(IO, dims=2)
    α = vadd[:, year-1959] ./ grossy[:, year-1959]
    β = ((I - diagm(1 .- α) * Ω)' * grossy[:, year-1959])
    β = max.(β,0.0) 
    β = β / sum(β)
    λ = (inv(I - diagm(1 .- α) * Ω)' * β)
    L = λ .* α

    return α, β, Ω, L, λ

end

data = loadInData();

grossy = reshape(data[:, 3], 46, 88)'
capital = reshape(data[:, 4], 46, 88)'; #NOMINAL CAPITAL
labor = reshape(data[:, 5], 46, 88)'; #NOMINAL LABOR
vadd = labor + capital; #NOMINAL VALUE;

removableSectors = [60, 80:88...];
grossy = grossy[setdiff(1:end, removableSectors), :];
capital = capital[setdiff(1:end, removableSectors), :];
labor = labor[setdiff(1:end, removableSectors), :];
vadd = vadd[setdiff(1:end, removableSectors), :];

grossSales = vec((sum(grossy, dims=2) .!= 0));
grossy = grossy[grossSales, :];
capital = capital[grossSales, :];
labor = labor[grossSales, :];
vadd = vadd[grossSales, :];

year = 1980

α, β, Ω, L, λ = getVariables(year)

size = 76

Ωlite = Ω[1:size, 1:size]
G = SimpleWeightedDiGraph(Ωlite)

d = Ω*(grossy[:,year - 1959]/sum(grossy[:,year - 1959]))

add_vertex!(G)

for i ∈ 1:(nv(G)-1)
    #TODO: REPLACE THE ONES WITH THE ACTUAL SALES
    add_edge!(G, i, size + 1, d[i])
end

adjacency_matrix(G)[:,77]

rem_edge!(G, size + 1, size + 1)

## ELASTICITIES
θ = 0.9

function changeProduction!(G, sector, amount, count)
    #println("Sector $sector Amount: $amount, Count $count")
    if count ≥ 10 || abs(amount - 1) ≤ 0.000001  #Eigenproduktion, can surely be estimated with elasticity of labor
        return
    end
    for nb in inneighbors(G, sector)
        w = get_weight(G, nb, sector)

        #println("Setting weight of Edge $nb -> $sector to $(w + w *(amount -1)), from $w")
        #println("Count: $count")
        change = (w * (amount - 1))
        add_edge!(G, nb, sector, w + change)
        #Again, there could be a elasticity included in the model
        changeProduction!(G, nb, 1 + change, count + 1)
    end
end


function shockConsumption!(G, sector::Integer, amount)
    w = get_weight(G, sector, nv(G))

    change = (w * (amount - 1))
    add_edge!(G, sector, nv(G), w + change)
    #println("Setting weight of edge 7 -> $(nv(G)) to $(amount * w)")
    changeProduction!(G, sector,1 + change, 0)
end

G2 = deepcopy(G)

function shock(sector,amount,G)
    G2 = deepcopy(G) 
    shockConsumption!(G2, sector, amount)
    x = sum(adjacency_matrix(G2)[1:76,1:76],dims=2)
    println(norm(adjacency_matrix(G) - adjacency_matrix(G2)))
    extrema(x)
end


y = LinRange(0,2,100)
f = [shock(55,x,G) for x in y]

f

shock(55,0,G)

edws = [get_weight(G2, e.src, e.dst) for e in edges(G2)]

graphplot(G2, edge_width=2 .* edws, arrow_size=5 .+ 5 .* edws)