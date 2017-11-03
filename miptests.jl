using JuMP
using Base.Test


# teste exemplo
# adicionado por Joaquim Garcia
function test1(solveMIP::Function)
    @testset "Teste examplo" begin
        m = Model()
        @variable(m, x >=0)
        @objective(m, :Min, x)

        sol = solveMIP(m)
        @test getobjectivevalue(m) == 0
        @test getvalue(x) == 0

        # TODO testar conteudo da struct "sol"
    end
end


#teste problema 1 da lista (mochila)
#adicionado por Guilherme Bodin
function test2(solveMIP::Function)
    @testset "Teste da Mochila" begin
        m = Model()
        @variable(m, x[i=1:3], Bin)
        @constraint(m, 6*x[1] + 5*x[2] + 5*x[3] <= 10)
        @objective(m, Max, 6*x[1] + 4*x[2] + 3*x[3])

        sol = solveMIP(m)
        @test getobjectivevalue(m) == 7
        @test getvalue(x) == [0, 1, 1]

        # TODO testar conteudo da struct "sol"
    end
end

#teste problema 6 da lista (Expansao da Producao)
#adicionado por Andrew Rosemberg
function test3(solveMIP::Function)
    Cinv = 13.16
    M = 200
    @testset "Teste da Expansao da Producao" begin
        model = Model()
        @variable(model, x[i=1:2]>=0)
        @variable(model, u, Bin)
        @objective(model, Max, 4*x[1] + 3*x[2] - u*Cinv)

        @constraint(model, 2*x[1] + 1*x[2] <= 4 +u*M)
        @constraint(model, 1*x[1] + 2*x[2] <= 4 +u*M)

        @constraint(model, 1*x[1] + 0.1*x[2] <= 4 +(1-u)*M)
        @constraint(model, 0.4*x[1] + 1*x[2] <= 4 +(1-u)*M)

        sol = solveMIP(model)
        @test getobjectivevalue(m) == 9.340000000000002
        @test getvalue(x) == [3.75, 2.5]
        @test getvalue(u) == 1

        # TODO testar conteudo da struct "sol"
    end
end


#teste problema da P1 (CVRP)
#adicionado por Eduardo Brito
function test4(solveMIP::Function)
    @testset "Teste CVRP" begin
        tam = 6
        carros = 2
        c = [ 0.162329   0.934386    0.27999    0.771633  0.93058    0.757016
        0.168044   0.516489    0.561025   0.793173  0.628186   0.893351
        0.450754   0.130342    0.0550682  0.330183  0.0140276  0.5666
        0.405403   0.00683533  0.828927   0.361477  0.265888   0.437654
        0.164714   0.00988625  0.470629   0.941743  0.343984   0.906805
        0.0105587  0.825022    0.540088   0.840939  0.137927   0.637206] #Custo na rota i,j
        b = [0.3286854210610859
        0.8747455782649833
        0.2997874087044681
        0.012584530553862994
        0.5477965855045668
        0.5671227617656462] #Demandas nos vertices
        b[1] = -sum(b[2:tam])
        k = [15.2168   9.21175  46.9288   2.46427  34.4648  25.4713
        84.4671  84.7527   84.7259  30.5085   55.8006  48.8416
        59.8025  87.1999   62.5169  57.5668   47.5803  66.2564
        97.9814  53.4096   72.0769   1.71473  41.5536  80.0004
        30.9586  59.9041   65.278   88.7064   55.0077  43.9505
        46.67    64.3936   32.2418  82.8831   38.0806  68.6481] #Capacidade na rota i,j

        m = Model()

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
         resp_x = [ 0.0  0.0       1.72233  0.579707  0.0      0.0
         0.0  0.0       0.0      0.0       0.0      0.0
         0.0  0.0       0.0      0.0       1.42254  0.0
         0.0  0.0       0.0      0.0       0.0      0.567123
         0.0  0.874746  0.0      0.0       0.0      0.0
         0.0  0.0       0.0      0.0       0.0      0.0]

         resp_y = [0.0  -0.0   1.0   1.0  -0.0   0.0
        1.0   0.0  -0.0  -0.0  -0.0  -0.0
        -0.0  -0.0   0.0   0.0   1.0  -0.0
        -0.0   0.0  -0.0   0.0  -0.0   1.0
        0.0   1.0  -0.0  -0.0   0.0  -0.0
        1.0  -0.0  -0.0  -0.0  -0.0   0.0]
        sol = solveMIP(m)
        @test getobjectivevalue(m) == 1.69179355
        @test getvalue(x) == resp_x
        @test getvalue(y) == resp_y

        # TODO testar conteudo da struct "sol"
    end
end


# teste Problema da Producao (PL)
# adicionado por Eduardo Brito
function test5(solveMIP::Function)
    @testset "Problema da Producao" begin
        m = Model(solver=GurobiSolver())
        @variable(m, x[1:2] >=0)
        @constraints(m, begin
        cons1, 2x[1] + x[2] <= 4
        cons2, x[1] + 2x[2] <= 4
        end)
        @objective(m, :Max, 4x[1] + 3x[2])

        sol = solve(m)
        @test getobjectivevalue(m) <= 9.33334
        @test getobjectivevalue(m) >= 9.33332
        @test getvalue(x) .<= [1.3333334;1.3333334]
        @test getvalue(x) .>= [1.3333332;1.3333332]

        # TODO testar conteudo da struct "sol"
    end
end
