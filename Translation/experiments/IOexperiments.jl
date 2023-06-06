using Graphs, GraphMakie, GLMakie, MAT, LinearAlgebra, SimpleWeightedGraphs, Printf

using GraphMakie.NetworkLayout

function loadInData()
    file = matopen("Translation/simulationData.mat")

    data = read(file, "data")

    close(file)

    return data
end

function getVariables(year)
    IO = data[data[:, 1].==year, :]
    temp = [8, 60, 62, 80:88...]
    IO = IO[setdiff(1:end, temp), setdiff(1:end, [1, 2, 3, 4, 5, 94])]
    IO = IO[:, setdiff(1:end, temp)]
    Ω = IO ./ sum(IO, dims=2)
    α = vadd[:, year-1959] ./ grossy[:, year-1959]
    β = ((I - diagm(1 .- α) * Ω)' * grossy[:, year-1959])
    @views β[β.<0] .= 0
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





α, β, Ω, L, λ = getVariables(1980)


Ωlite = Ω[1:10, 1:10]
G = SimpleWeightedDiGraph(Ωlite)


add_vertex!(G)

for i ∈ 1:(nv(G)-1)
    add_edge!(G, i, 11, 1)
end

rem_edge!(G,11,11)


function increaseProduction!(G, sector, amount)
    if abs(amount - 1)  ≤ 0.05
        return
    end
    for nb in inneighbors(G, sector)
        w = get_weight(G, nb, sector)

        println("Setting weight of Edge $nb -> $sector to $(w * amount), from $w")
        add_edge!(G, nb, sector, w * amount)

        increaseProduction!(G, nb, 1 + w * (amount - 1))
    end
end

function shockConsumption!(G, sector::Integer, amount)
    w = (get_weight(G, sector, 5))

    add_edge!(G, sector, 5, amount * w)

    increaseProduction!(G, sector, w * amount)
end

G2 = deepcopy(G)


shockConsumption!(G2, 9, 1.3)

adjacency_matrix(G) - adjacency_matrix(G2)

edws = [get_weight(G2,e.src,e.dst) for e in edges(G2)]

string_x = [@sprintf("%5.3f",x) for x in edws]

graphplot(G2,edge_width =2 .* edws, arrow_size = 5 .+ 5 .* edws,elabels = string_x)




