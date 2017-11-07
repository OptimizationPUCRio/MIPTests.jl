using JuMP
using MathProgBase
using Base.Test

#≈

# teste exemplo
# adicionado por Joaquim Garcia
function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste examplo" begin
        m = Model(solver = solver)
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
function test2(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste da Mochila" begin
        m = Model(solver = solver)
        @variable(m, x[i=1:3], Bin)
        @constraint(m, 6*x[1] + 5*x[2] + 5*x[3] <= 10)
        @objective(m, Max, 6*x[1] + 4*x[2] + 3*x[3])

        sol = solveMIP(m)
        @test getobjectivevalue(m) == 7
        @test getvalue(x) == [0, 1, 1]

        # TODO testar conteudo da struct "sol"
    end
end

# teste Sudoku
# adicionado por Raphael Saavedra
function testSudoku(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste Sudoku" begin
        n = 9
        model = Model()
        @variable(model, x[i in 1:n, j in 1:n, k in 1:n], Bin)

        fixas = [(1,3,4), (1,5,6), (1,9,2), (2,1,8), (2,3,5), (2,6,2), (2,8,3),
                (3,5,3), (3,8,6), (4,2,2), (4,3,8), (5,6,4), (6,1,7), (6,5,5),
                (6,9,9), (7,3,2), (7,6,1), (8,2,7), (8,5,4), (8,7,9), (8,8,5),
                (9,1,6), (9,8,4)]
        for idx in fixas
            @constraint(model, x[idx...] == 1)
        end
        @constraint(model, [j in 1:n, k in 1:n], sum(x[:,j,k]) == 1)
        @constraint(model, [i in 1:n, k in 1:n], sum(x[i,:,k]) == 1)
        @constraint(model, [i in 1:n, j in 1:n], sum(x[i,j,:]) == 1)
        @constraint(model, [p in [0,3,6], q in [0,3,6], k in 1:n], sum(sum(x[i+p,j+q,k] for i in 1:3) for j in 1:3) == 1)
        @objective(model, Min, 0)

        sol = solveMIP(model)
        @test getobjectivevalue(m) == 0
        @test sum(getvalue(x)) == 81
        @test getvalue(x[1,1,1]) == 0
        @test getvalue(x[1,1,3]) == 1
        @test getvalue(x[8,1,1]) == 0
        @test getvalue(x[8,1,9]) == 1
        @test sum(getvalue(x[6,6,4:9])) == 0
        @test getvalue(x[6,6,3]) == 1
    end
end

#teste problema 6 da lista (Expansao da Producao)
#adicionado por Andrew Rosemberg
function test3(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    Cinv = 13.16
    M = 200
    @testset "Teste da Expansao da Producao" begin
        model = Model(solver = solver)
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

=======
# teste mochila binária infeasible
# adicionado por Raphael Saavedra
function testInfeasibleKnapsack(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste Mochila infeasible" begin
        m = Model(solver = solver)
        @variable(m, x[i=1:3], Bin)
        @constraint(m, 6*x[1] + 5*x[2] + 5*x[3] <= 5)
        @constraint(m, x[1] == 1)
        @objective(m, Max, 6*x[1] + 4*x[2] + 3*x[3])

        sol = solveMIP(m)
        @test sol.ext[:status] == :Infeasible

        >>>>>>> origin/Branch-Raphael
    end
end


#teste problema da P1 (CVRP)
#adicionado por Eduardo Brito
function test_P1_Brito(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
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

        m = Model(solver=solver)

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
        @test getobjectivevalue(m) ≈ 1.69179355 atol=exp10(-5)
        @test getvalue(x) == resp_x
        @test getvalue(y) == resp_y

        # TODO testar conteudo da struct "sol"
    end
end


# teste Problema da Producao (PL)
# adicionado por Eduardo Brito
function test_PL_Simples_Brito(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Problema da Producao" begin
        m = Model(solver=solver())
        @variable(m, x[1:2] >=0)
        @constraints(m, begin
        cons1, 2x[1] + x[2] <= 4
        cons2, x[1] + 2x[2] <= 4
        end)
        @objective(m, :Max, 4x[1] + 3x[2])

        sol = solveMIP(m)
        @test getobjectivevalue(m) ≈ 9.33334 atol = exp10(-5)
        @test getvalue(x) ≈ [1.3333334;1.3333334] atol = exp10(-5)

        # TODO testar conteudo da struct "sol"
    end
end


# teste Pl Infeasible
# adicionado por Eduardo Brito
function test_PL_Infeasible_Brito(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Problema Infeasible" begin
        m = Model(solver=solver)
        @variable(m, x[1:2] >=0)
        @constraints(m, begin
        cons1, 2x[1] + x[2] <= 4
        cons2, x[1] + 2x[2] <= 4
        cons3, x[1] + x[2] >= 5
        end)
        @objective(m, :Max, 4x[1] + 3x[2])

        sol = solveMIP(m)

        @test sol.ext[:status] == :Infeasible

        # TODO testar conteudo da struct "sol"
    end
end


# teste Pl Unbounded
# adicionado por Eduardo Brito
function test_PL_Unbounded_Brito(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Problema Unbounded" begin
        m = Model(solver=solver)
        @variable(m, x[1:2] >=0)
        @objective(m, :Max, 4x[1] + 3x[2])

        sol = solveMIP(m)
        @test sol.ext[:status] == :Unbounded

        # TODO testar conteudo da struct "sol"
    end
end


# teste MIP (minimal ~5 binarias)
# adicionado por Eduardo Brito
function test_MIP_Minimal_Brito(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "MIP Minimal" begin
        m = Model(solver=solver)
        @variable(m, x[1:5] >=0, Bin)
        @variable(m, y[1:5] >= 0)
        @constraints(m, begin
        cons1, sum(x) <= 4.5
        cons2, y[1] <= 10(x[1])
        cons3, y[2] <= 10(x[2])
        cons4, y[3] <= 10(x[3])
        cons5, y[4] <= 10(x[4])
        cons6, y[5] <= 10(x[5])
        end)
        @objective(m, :Max, 5y[1] + 4y[2] + 3y[3] + 2y[4] + 1y[5])

        sol = solveMIP(m)
        @test getobjectivevalue(m) == 140
        @test getvalue(x) == [1;1;1;1;0]
        @test getvalue(y) ≈ [10;10;10;10;0] atol = exp10(-5)

        # TODO testar conteudo da struct "sol"
    end
end

# teste MIP Pequeno (~50 binarias) ~ The Assignment problem:
# adicionado por Eduardo Brito
function test_MIP_Pequeno_Brito(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste MIP Pequeno, Assignment problem" begin
    n  = 8
    if true
        c = [0.445321 0.499462 0.409753 0.471825 0.1172 0.820595 0.629809 0.333445; 0.197025 0.160481 0.00865311 0.355901 0.137367 0.199186 0.718575 0.716486;
         0.654497 0.904598 0.321483 0.171736 0.050554 0.254487 0.540093 0.724331; 0.254369 0.593379 0.205166 0.288702 0.499699 0.308233 0.869406 0.353904;
         0.854515 0.00978121 0.520072 0.985762 0.72076 0.317384 0.268573 0.315585; 0.0212753 0.754076 0.753672 0.158407 0.212617 0.403343 0.71157 0.17261;
         0.651835 0.24596 0.700141 0.989018 0.723494 0.236829 0.891181 0.568245; 0.257637 0.883802 0.0252095 0.0273074 0.450492 0.560833 0.820861 0.893546]
    end
    m = Model(solver = solver)
    @variable(m, x[i=1:n,j=1:n] >=0, Bin)
    @objective(m, :Min, sum(c[i,j]*x[i,j] for i = 1:n, j = 1:n))
    @constraints(m,begin
    linhas[i=1:n], sum(x[i,j] for j = 1:n) == 1
    colunas[j=1:n], sum(x[i,j] for i = 1:n) == 1
    end)

    sol = solveMIP(m)
    @test getobjectivevalue(m) ≈ 1.264 atol = exp10(-5)
    @test abs.(getvalue(x)) == abs.([0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0;
                                    0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0;
                                    0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0;
                                    0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0;
                                    0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0;
                                    1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                                    0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0;
                                    0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0])


        # TODO testar conteudo da struct "sol"
    end
end
