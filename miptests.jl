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
        model = Model(solver = solver)
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

# teste n-K robusto - trabalho da P1 - caso com 10 geradores e K = 2
# adicionado por Raphael Saavedra
function testRobustCCUC(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
  #------------------------------------------------------------------------------
  # Parameters
  T = collect(1:24) # periods
  N = collect(1:10) # generators
  p0 = zeros(10) # initial power output
  v0 = zeros(10) # initial on/off state
  D = 0.5*[700; 750; 850; 950; 1000; 1100; 1150; 1200; 1300; 1400; 1450; 1500; 1400;
            1300; 1200; 1050; 1000; 1100; 1200; 1400; 1300; 1100; 900; 800] # demand
  K = 2 # n-K security criterion
  Cf = [1000; 970; 700; 680; 450; 370; 480; 660; 665; 670] # fixed cost
  Cl = [16.19; 17.26; 16.6; 16.5; 19.7; 22.26; 27.74; 25.92; 27.27; 27.79] # linear cost
  Cs = 0.08*Cl # spinning reserve cost
  Cns = 0.1*Cl # non-spinning reserve cost
  Pmax = [455; 455; 130; 130; 162; 80; 85; 55; 55; 55] # generator capacity
  Pmin = [150; 150; 20; 20; 25; 20; 25; 10; 10; 10] # minimum power output
  RSmax = Pmax # maximum spinning reserve
  RNSmax = Pmax # maximum non-spinning reserve
  RD = Pmax # ramp-down limit
  RU = Pmax # ramp-up limit
  SD = RD # shutdown ramp limit
  SU = RU # startup ramp limit
  #------------------------------------------------------------------------------
  @testset "Teste Robust n-K Unit Commitment" begin
    # Model formulation
    m = Model(solver = solver)
    #----------------------------------------------------------------------------
    # Variables
    @variable(m, p[1:N[end], 1:T[end]] >= 0) # power output
    @variable(m, v[1:N[end], 1:T[end]], Bin) # 1 if generator is on, 0 otherwise
    @variable(m, 0 <= rs[i = 1:N[end], 1:T[end]] <= RSmax[i]) # spinning reserve
    @variable(m, rns[1:N[end], 1:T[end]] >= 0) # non-spinning reserve
    @variable(m, vns[1:N[end], 1:T[end]], Bin) # 1 if generator provides non-spinning reserve, 0 otherwise
    @variable(m, y[1:T[end]] >= 0) # dual variable of the n-K constraint
    @variable(m, z[1:N[end], 1:T[end]] >= 0) # dual variable of the upper bound
    #----------------------------------------------------------------------------
    # Constraints
    @constraint(m, [t in T], sum(p[i,t] for i in N) == D[t])
    @constraint(m, [t in T, i in N], p[i,t] >= Pmin[i]*v[i,t])
    @constraint(m, [t in T, i in N], p[i,t] <= Pmax[i]*v[i,t])
    @constraint(m, [t in T, i in N], p[i,t] + rs[i,t] <= Pmax[i]*v[i,t])
    @constraint(m, [t in T, i in N], rns[i,t] >= Pmin[i]*vns[i,t])
    @constraint(m, [t in T, i in N], rns[i,t] <= RNSmax[i]*vns[i,t])
    @constraint(m, [t in T, i in N], v[i,t] + vns[i,t] <= 1)
    @constraint(m, [t in 2:T[end], i in N], p[i,t-1] <= p[i,t] + RD[i]*v[i,t] + SD[i]*(v[i,t-1]-v[i,t]) + Pmax[i]*(1-v[i,t-1]))
    @constraint(m, [t in 2:T[end], i in N], p[i,t] <= p[i,t-1] + RU[i]*v[i,t-1] + SU[i]*(v[i,t]-v[i,t-1]) + Pmax[i]*(1-v[i,t]))
    @constraint(m, [t in 2:T[end], i in N], p[i,t] + rs[i,t] <= p[i,t-1] + RU[i]*v[i,t-1] + SU[i]*(v[i,t]-v[i,t-1]) + Pmax[i]*(1-v[i,t]))
    @constraint(m, [t in 2:T[end], i in N], rns[i,t] <= p[i,t-1] + RU[i]*v[i,t-1] + RNSmax[i]*(1-v[i,t-1]))
    @constraint(m, [t in 2:T[end], i in N], rns[i,t] <= SU[i]*(vns[i,t]-v[i,t-1]) + RNSmax[i]*(1-(vns[i,t]-v[i,t-1])))
    @constraint(m, [t in T], (N[end]-K)*y[t] - sum(z[i,t] for i in N) >= D[t])
    @constraint(m, [t in T, i in N], y[t] - z[i,t] <= p[i,t] + rs[i,t] + rns[i,t])

    # Constraints regarding initial state
    @constraint(m, [t in 1, i in N], p0[i] <= p[i,t] + RD[i]*v[i,t] + SD[i]*(v0[i]-v[i,t]) + Pmax[i]*(1-v0[i]))
    @constraint(m, [t in 1, i in N], p[i,t] <= p0[i] + RU[i]*v0[i] + SU[i]*(v[i,t]-v0[i]) + Pmax[i]*(1-v[i,t]))
    @constraint(m, [t in 1, i in N], p[i,t] + rs[i,t] <= p0[i] + RU[i]*v0[i] + SU[i]*(v[i,t]-v0[i]) + Pmax[i]*(1-v[i,t]))
    @constraint(m, [t in 1, i in N], rns[i,t] <= p0[i] + RU[i]*v0[i] + RNSmax[i]*(1-v0[i]))
    @constraint(m, [t in 1, i in N], rns[i,t] <= SU[i]*(vns[i,t]-v0[i]) + RNSmax[i]*(1-(vns[i,t]-v0[i])))
    #----------------------------------------------------------------------------
    # Objective function
    @objective(m, Min, sum(Cf[i]*v[i,t] for i in N, t in T) + sum(Cl[i]*p[i,t] for i in N, t in T) +
                            sum(Cs[i]*rs[i,t] for i in N, t in T) + sum(Cns[i]*rns[i,t] for i in N, t in T))
  #------------------------------------------------------------------------------
    sol = solveMIP(m)

    @test getobjectivevalue(m) ≈ 289892.9539 atol=1e-3
    @test getvalue(p[:,24]) == [400; zeros(9)]
    @test getvalue(v[1,:]) == ones(24)
    @test getvalue(v[:,1]) == [1; zeros(9)]
    @test getvalue(v[:,4]) == [1; zeros(4); 1; zeros(4)]
    @test sum(getvalue(v)) == 43

  end
end

# teste Caminho mais curto
# adicionado por Carlos
function testCaminho(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste caminho mais curto" begin

        m = Model(solver = solver)
        @variable(m, f[i in 1:6, j in 1:6], Bin)

        A = [0 1 1 0 0 0
             0 0 0 1 0 0
             0 1 0 1 1 0
             0 0 0 0 1 1
             0 0 0 0 0 1
             0 0 0 0 0 0]

        c = [0 2 2 0 0 0
             0 0 0 3 0 0
             0 1 0 1 3 0
             0 0 0 0 1 1
             0 0 0 0 0 2
             0 0 0 0 0 0]

        b = [1;0;0;0;0;-1]

        @constraint(m,[v=1:6], sum(A[v,j]*x[v,j] for j=1:6) - sum(A[i,v]*x[i,v] for i=1:6) == b[v])

        @objective(m, Min, sum(A[i,j]*c[i,j]*x[i,j] for i=1:6, j=1:6))

        sol = solveMIP(m)
        @test getobjectivevalue(m) == 4
        @test getvalue(x[1,3]) == 1
        @test getvalue(x[3,4]) == 1
        @test getvalue(x[4,6]) == 1
        @test sum(getvalue(x)) == 3
    end
end

