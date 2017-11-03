using JuMP, Gurobi, LightGraphs, GraphPlot, Gadfly

m = Model(solver=GurobiSolver())

tam = 20
carros = 3
c = rand(tam,tam) #Custo na rota i,j
b = rand(tam) #Demandas nos vertices
b[1] = -sum(b[2:tam])
k = 100*rand(tam,tam) #Capacidade na rota i,j


@variable(m,x[1:tam,1:tam] >= 0)
@variable(m,y[1:tam,1:tam], Bin)

@constraints(m,begin
saida[i=2:tam],sum(y[i,j] for j = 1:tam if j!=i) == 1
chegada[j=2:tam],sum(y[i,j] for i = 1:tam if j!=i) == 1
s_depot, sum(y[i,1] for i = 2:tam) == carros
c_depot, sum(y[1,j] for j = 2:tam) == carros
end)

@constraints(m,begin
capacidade_rota[i=1:tam,j=1:tam], x[i,j] <= k[i,j]*y[i,j]
capacidade_carro[i=1:tam], x[1,i] <= sum(b[2:tam])/carros + sqrt(var(b))
demanda[i=1:tam], sum(x[j,i] for j =1:tam if j!= i) - sum(x[i,j] for j=1:tam if j!=i) == b[i]
ciclos, sum(y[i,i] for i =1:tam) == 0
end)

@objective(m, Min, sum(sum(c[i,j]*y[i,j] for i = 1:tam) for j = 1:tam))

solve(m)

getvalue(x)
g = DiGraph(abs.(getvalue(y)))
membership = 2*ones(Int8,nv(g));
membership[1] = 1;
nodecolor = [colorant"orange",colorant"lightseagreen"];
nodefillc = nodecolor[membership];
nodelabel = collect(1:1:nv(g));
gplot(g,nodelabel = nodelabel,nodefillc=nodefillc)
ind = simplecycles_hadwick_james(g) #Achando os ciclos
