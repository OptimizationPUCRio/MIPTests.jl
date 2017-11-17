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


# adicionado por Bianca Lacê
function test3_MIP_minimal_Bianca(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    model = Model(solver = solver)
    c=[1 ; 1 ; 3 ; 3 ; 2]
    M=15
    k=3
    testresult = @testset "Teste MIP minimo" begin
        n = length(c)
        @variables(model, begin
          z[j=1:n], Bin
          y[i=1:n] >= 0
        end)
        @constraint(model,  constrain[i=1:n], y[i]>= z[i]*M)
        @constraint(model,  sum(z[j] for j=1:n) == k)
        @objective(model, Min, sum(c[i]*y[i] for i=1:n))

        solveMIP(model)

        @test getobjectivevalue(model) == 60
        @test getvalue(z) == [1.00;1.00;0.00;0.00;1.00]
        @test getvalue(y) == [15.0;15.0;0.00;0.00;15.0]

    end
    setoutputs!(model,solution,testresult)
    return solution
end
