###testes
include("C:/Users/Bianca/Desktop/otimização/Prog Inteira/github/branchnbound.jl")

using JuMP
using MathProgBase
using Base.Test

type MIPSolution
    name::String
    pass::Bool
    objective::Float64
    bestbound::Float64
    time::Float64
    nodes::Int
    intsols::Int
    status::Symbol
    function MIPSolution()
        new("",false,NaN,NaN,Inf,-1,-1,:unsolved)
    end
end

function setoutputs!(m,sol::MIPSolution, test)
    sol.pass = !test.anynonpass
    sol.objective = getobjectivevalue(m)
    sol.bestbound = m.objBound
    if haskey(m.ext,:time)
        sol.time = m.ext[:time]
    end
    if haskey(m.ext,:nodes)
        sol.nodes = m.ext[:nodes]
    end
    if haskey(m.ext,:intsols)
        sol.intsols = m.ext[:intsols]
    end
    if haskey(m.ext,:status)
        sol.status = m.ext[:status]
    end
    return nothing
end


#≈

# teste problema de atribuição de tarefas
# adicionado por Bianca Lacê
function test1_Bianca(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    model = Model(solver = solver)
    C=[2 1 3; 1 2 5; 3 2 1]
    testresult = @testset "Teste atribuição de tarefas" begin
      l,n=size(C)
      @variable(m, X[1:l,1:n], Bin)
      @constraints(m, begin
        constrain[i=1:l], sum(X[i,j] for j=1:n) == 1
        constrain[j=1:n], sum(X[i,j] for i=1:l) == 1
        sum(sum(C[i,j]*X[i,j] for j=1:n) for i=1:l) <= 9
      end)
      @objective(m, :Min, sum(sum(C[i,j]*X[i,j] for j=1:n) for i=1:l))

      solveMIP(m)
      @test getobjectivevalue(m) == 3
      @test getvalue(X) == [0 1 0; 1 0 0; 0 0 1]

    end
    setoutputs!(m,solution,testresult)
    return solution
end

# teste problema de atribuição de tarefas inviavel
# adicionado por Bianca Lacê
function test2_Bianca(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    C=[10 20 30; 15 12 20; 10 25 13]
    testresult = @testset "Teste atribuição de tarefas" begin
        l,n=size(C)
        @variable(m, X[1:l,1:n], Bin)
        @constraints(m, begin
          constrain[i=1:l], sum(X[i,j] for j=1:n) == 1
          constrain[j=1:n], sum(X[i,j] for i=1:l) == 1
          sum(sum(C[i,j]*X[i,j] for j=1:n) for i=1:l) <= 9
        end)
        @objective(m, :Min, sum(sum(C[i,j]*X[i,j] for j=1:n) for i=1:l))

        solveMIP(m)
        @test m.ext[:status] == :Infeasible

    end
    setoutputs!(m,solution,testresult)
    return solution
end
