using JuMP
using MathProgBase
using Base.Test


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
        @test getobjectivevalue(model) ≈ 9.340000000000002 atol=1E-07
        @test getvalue(x) ≈ [3.75, 2.5] atol=1E-07
        @test getvalue(u) ≈ 1 atol=1E-07

        # TODO testar conteudo da struct "sol"
    end
end

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
        @test sol == :Infeasible
        @test m.ext[:status] == :Infeasible
    end
end
