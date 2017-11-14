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
        m = Model(solver = solver)
        @variable(m, x[i in 1:n, j in 1:n, k in 1:n], Bin)

        fixas = [(1,3,4), (1,5,6), (1,9,2), (2,1,8), (2,3,5), (2,6,2), (2,8,3),
                (3,5,3), (3,8,6), (4,2,2), (4,3,8), (5,6,4), (6,1,7), (6,5,5),
                (6,9,9), (7,3,2), (7,6,1), (8,2,7), (8,5,4), (8,7,9), (8,8,5),
                (9,1,6), (9,8,4)]
        for idx in fixas
            @constraint(m, x[idx...] == 1)
        end
        @constraint(m, [j in 1:n, k in 1:n], sum(x[:,j,k]) == 1)
        @constraint(m, [i in 1:n, k in 1:n], sum(x[i,:,k]) == 1)
        @constraint(m, [i in 1:n, j in 1:n], sum(x[i,j,:]) == 1)
        @constraint(m, [p in [0,3,6], q in [0,3,6], k in 1:n], sum(sum(x[i+p,j+q,k] for i in 1:3) for j in 1:3) == 1)
        @objective(m, Min, 0)

        sol = solveMIP(m)
        @test getobjectivevalue(m) == 0
        @test sum(getvalue(x)) == 81
        @test getvalue(x[1,1,1]) == 0
        @test getvalue(x[1,1,3]) == 1
        @test getvalue(x[8,1,1]) == 1
        @test getvalue(x[8,1,9]) == 0
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

        @test m.ext[:status] == :Infeasible


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

        solveMIP(m)
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
        m = Model(solver=solver)
        @variable(m, x[1:2] >=0)
        @constraints(m, begin
        cons1, 2x[1] + x[2] <= 4
        cons2, x[1] + 2x[2] <= 4
        end)
        @objective(m, :Max, 4x[1] + 3x[2])

        solveMIP(m)
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

        solveMIP(m)

        @test m.ext[:status] == :Infeasible

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

        solveMIP(m)
        @test m.ext[:status] == :Unbounded

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

        solveMIP(m)
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

    solveMIP(m)
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

        @test m.ext[:status] == :Optimal

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
        @variable(m, x[i in 1:6, j in 1:6], Bin)

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


#teste problema 6 da lista modificado 1 (Expansao da Producao Unbounded)
#adicionado por Andrew Rosemberg
function test3_2(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    Cinv = 13.16
    M = 200
    @testset "Teste da Expansao da Producao Unbounded" begin
        model = Model(solver = solver)
        @variable(model, x[i=1:2]>=0)
        @variable(model, u, Bin)
        @objective(model, Max, 4*x[1] + 3*x[2] - u*Cinv)

        @constraint(model, 2*x[1] + 1*x[2] >= 4 +u*M)
        @constraint(model, 1*x[1] + 2*x[2] >= 4 +u*M)

        @constraint(model, 1*x[1] + 0.1*x[2] >= 4 +(1-u)*M)
        @constraint(model, 0.4*x[1] + 1*x[2] >= 4 +(1-u)*M)

        solveMIP(model)

        @test model.ext[:status] == :Unbounded
    end
end


#teste problema 6 da lista modificado 2 (Expansao da Producao Infeasible)
#adicionado por Andrew Rosemberg
function test3_3(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    Cinv = 13.16
    M = 200
    @testset "Teste da Expansao da Infeasible" begin
        model = Model(solver = solver)
        @variable(model, x[i=1:2]>=0)
        @variable(model, u, Bin)
        @objective(model, Max, 4*x[1] + 3*x[2] - u*Cinv)

        @constraint(model, 2*x[1] + 1*x[2] <= 4 -u*M)
        @constraint(model, 1*x[1] + 2*x[2] <= 4 -u*M)

        @constraint(model, 1*x[1] + 0.1*x[2] <= 4 -(1-u)*M)
        @constraint(model, 0.4*x[1] + 1*x[2] <= 4 -(1-u)*M)

        solveMIP(model)

        @test model.ext[:status] == :Infeasible
    end
end

#teste Feature selection pequeno (Viavel)
#adicionado por Andrew Rosemberg
function test_feature_selection_pequeno_viavel(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    srand(2)
    numpossiblevar = 50
    numvar = 40
    numconstraints = 50
    datamatrix = rand(numconstraints,numpossiblevar)
    weights = [rand(numvar);zeros(numpossiblevar-numvar)]
    index_vector = datamatrix*weights
    maxweight  = maximum(weights)+1
    @testset "Teste Feature selection pequeno Viavel" begin
        model = Model(solver = solver)
        @objective(model, Max, 0)

        @variable(model, w[i=1:numpossiblevar]>=0)
        @variable(model, u[i=1:numpossiblevar], Bin)

        @constraint(model, sum(u) == numvar)

        @constraints(model,begin
          dummy[i=1:numpossiblevar], w[i] <= u[i]*maxweight
        end)

        @constraints(model,begin
          linhas[i=1:numconstraints], (datamatrix*w)[i] == index_vector[i]
        end)

        solveMIP(model)

        utest = [ones(numvar);zeros(numpossiblevar-numvar)]
        @test utest == getvalue(u)
        @test getvalue(w) ≈ weights atol=1E-07
    end
end

#teste Feature selection medio (Viavel)
#adicionado por Andrew Rosemberg
function test_feature_selection_medio(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    srand(1000)
    numpossiblevar = 500
    numvar = 15
    numconstraints = 100
    datamatrix = rand(numconstraints,numpossiblevar)
    weights = [rand(numvar);zeros(numpossiblevar-numvar)]
    index_vector = datamatrix*weights
    maxweight  = maximum(weights)+1
    @testset "Teste Feature selection pequeno Viavel" begin
        model = Model(solver = solver)
        @objective(model, Max, 0)

        @variable(model, w[i=1:numpossiblevar]>=0)
        @variable(model, u[i=1:numpossiblevar], Bin)

        @constraint(model, sum(u) == numvar)

        @constraints(model,begin
          dummy[i=1:numpossiblevar], w[i] <= u[i]*maxweight
        end)

        @constraints(model,begin
          linhas[i=1:numconstraints], (datamatrix*w)[i] == index_vector[i]
        end)

        solveMIP(model)

        utest = [ones(numvar);zeros(numpossiblevar-numvar)]
        @test utest == getvalue(u)
        @test getvalue(w) ≈ weights atol=1E-07
    end
end

#teste Feature selection grande (Viavel)
#adicionado por Andrew Rosemberg
function test_feature_selection_grande(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    srand(1000)
    numpossiblevar = 5000
    numvar = 15
    numconstraints = 100
    datamatrix = rand(numconstraints,numpossiblevar)
    weights = [rand(numvar);zeros(numpossiblevar-numvar)]
    index_vector = datamatrix*weights
    maxweight  = maximum(weights)+1
    @testset "Teste Feature selection pequeno Viavel" begin
        model = Model(solver = solver)
        @objective(model, Max, 0)

        @variable(model, w[i=1:numpossiblevar]>=0)
        @variable(model, u[i=1:numpossiblevar], Bin)

        @constraint(model, sum(u) == numvar)

        @constraints(model,begin
          dummy[i=1:numpossiblevar], w[i] <= u[i]*maxweight
        end)

        @constraints(model,begin
          linhas[i=1:numconstraints], (datamatrix*w)[i] == index_vector[i]
        end)

        solveMIP(model)

        utest = [ones(numvar);zeros(numpossiblevar-numvar)]
        @test utest == getvalue(u)
        @test getvalue(w) ≈ weights atol=1E-07
    end
end

#teste Feature selection pequeno (Inviavel)
#adicionado por Andrew Rosemberg
function test_feature_selection_pequeno_inviavel(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    srand(2)
    numpossiblevar = 50
    numvar = 40
    numconstraints = 50
    datamatrix = rand(numconstraints,numpossiblevar)
    weights = [rand(numvar);zeros(numpossiblevar-numvar)]
    index_vector = -datamatrix*weights
    maxweight  = maximum(weights)+1
    @testset "Teste Feature selection pequeno Viavel" begin
        model = Model(solver = solver)
        @objective(model, Max, 0)

        @variable(model, w[i=1:numpossiblevar]>=0)
        @variable(model, u[i=1:numpossiblevar], Bin)

        @constraint(model, sum(u) == numvar)

        @constraints(model,begin
          dummy[i=1:numpossiblevar], w[i] <= u[i]*maxweight
        end)

        @constraints(model,begin
          linhas[i=1:numconstraints], (datamatrix*w)[i] == index_vector[i]
        end)

        solveMIP(model)

        @test model.ext[:status] == :Infeasible
    end
end

#teste problema ucla : https://www.math.ucla.edu/~tom/LP.pdf pg 9 (PL unbounded)
#adicionado por Andrew Rosemberg
function teste_PL_andrew_unbounded(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste PL viavel da ucla" begin
        model = Model(solver = solver)
        @variable(model, x[i=1:3]>=0)
        @variable(model, x4)
        @objective(model, Max, x[1] + 3*x[2] + 4*x[3] + 2*x4 +5)

        @constraint(model, 4*x[1] + 2*x[2] +x[3] + 3*x4 <= 10)
        @constraint(model, x[1] - x[2] + 2*x[3] == 2)
        @constraint(model, x[1] + x[2] + x[3] + x4 >= 1)

        solveMIP(model)

        @test model.ext[:status] == :Unbounded
    end
end

#teste problema ucla modificado (PL viavel)
#adicionado por Andrew Rosemberg
function teste_PL_andrew_viavel(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste PL viavel da ucla" begin
        model = Model(solver = solver)
        @variable(model, x[i=1:3]>=0)
        @variable(model, x4>=0)
        @objective(model, Max, x[1] + 3*x[2] + 4*x[3] + 2*x4 +5)

        @constraint(model, 4*x[1] + 2*x[2] +x[3] + 3*x4 <= 10)
        @constraint(model, x[1] - x[2] + 2*x[3] == 2)
        @constraint(model, x[1] + x[2] + x[3] + x4 >= 1)

        solveMIP(model)

        @test getobjectivevalue(model) == 27
        @test getvalue(x) ≈ [0.0;3.6;2.8] atol=1E-07
        @test getvalue(x4) == 0
    end
end

#teste problema ucla modificado 2 (PL inviavel)
#adicionado por Andrew Rosemberg
function teste_PL_andrew_inviavel(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste PL viavel da ucla" begin
        model = Model(solver = solver)
        @variable(model, x[i=1:3]>=0)
        @variable(model, x4>=0)
        @objective(model, Max, x[1] + 3*x[2] + 4*x[3] + 2*x4 +5)

        @constraint(model, 4*x[1] + 2*x[2] +x[3] + 3*x4 <= -10)
        @constraint(model, x[1] - x[2] + 2*x[3] == 2)
        @constraint(model, x[1] + x[2] + x[3] + x4 >= 1)

        solveMIP(model)

        @test model.ext[:status] == :Infeasible
    end
end
# teste MIP unbounded
# adicionado por Raphael Saavedra
function testUnboundedKnapsack(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste Mochila unbounded" begin
        m = Model(solver = solver)
        @variable(m, x[i=1:3], Bin)
        @variable(m, y >= 0)
        @constraint(m, 6*x[1] + 5*x[2] + 5*x[3] <= 5)
        @objective(m, Max, 6*x[1] + 4*x[2] + 3*x[3] + y)

        sol = solveMIP(m)

        @test m.ext[:status] == :Unbounded
    end
end

# teste Infeasible Unit Commitment
# adicionado por Raphael Saavedra
function testInfeasibleUC(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
  #------------------------------------------------------------------------------
  # Parameters
  T = collect(1:10) # periods
  N = collect(1:5) # generators
  p0 = zeros(N[end]) # initial power output
  v0 = zeros(N[end]) # initial on/off state
  D = [400; 390; 380; 370; 360; 350; 340; 330; 320; 500] # demand
  Cf = [100; 100; 100; 100; 100] # fixed cost
  Cl = [10; 20; 30; 40; 50] # linear cost
  Pmax = [100; 100; 100; 100; 100] # generator capacity
  Pmin = [10; 10; 10; 10; 10] # minimum power output
  RD = [10; 20; 30; 40; 50] # ramp-down limit
  RU = [10; 25; 30; 40; 50] # ramp-up limit
  SD = RD # shutdown ramp limit
  SU = RU # startup ramp limit
  #------------------------------------------------------------------------------
  @testset "Teste Infeasible Unit Commitment" begin

    m = Model(solver = solver)
    #----------------------------------------------------------------------------
    # Variables
    @variable(m, p[1:N[end], 1:T[end]] >= 0) # power output
    @variable(m, v[1:N[end], 1:T[end]], Bin) # 1 if generator is on, 0 otherwise
    #----------------------------------------------------------------------------
    # Constraints
    @constraint(m, [t in T], sum(p[i,t] for i in N) == D[t])
    @constraint(m, [t in T, i in N], p[i,t] >= Pmin[i]*v[i,t])
    @constraint(m, [t in T, i in N], p[i,t] <= Pmax[i]*v[i,t])
    @constraint(m, [t in 2:T[end], i in N], p[i,t-1] <= p[i,t] + RD[i]*v[i,t] + SD[i]*(v[i,t-1]-v[i,t]) + Pmax[i]*(1-v[i,t-1]))
    @constraint(m, [t in 2:T[end], i in N], p[i,t] <= p[i,t-1] + RU[i]*v[i,t-1] + SU[i]*(v[i,t]-v[i,t-1]) + Pmax[i]*(1-v[i,t]))
    # Constraints regarding initial state
    @constraint(m, [t in 1, i in N], p0[i] <= p[i,t] + RD[i]*v[i,t] + SD[i]*(v0[i]-v[i,t]) + Pmax[i]*(1-v0[i]))
    @constraint(m, [t in 1, i in N], p[i,t] <= p0[i] + RU[i]*v0[i] + SU[i]*(v[i,t]-v0[i]) + Pmax[i]*(1-v[i,t]))
    #----------------------------------------------------------------------------
    # Objective function
    @objective(m, Min, sum(Cf[i]*v[i,t] for i in N, t in T) + sum(Cl[i]*p[i,t] for i in N, t in T))
  #------------------------------------------------------------------------------
    sol = solveMIP(m)

    @test m.ext[:status] == :Infeasible

  end
end

# teste PL simples
# adicionado por Raphael Saavedra
function test_PL_Simples_Raphael(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste PL simples" begin
        m = Model(solver = solver)
        @variable(m, x[i=1:3] >= 0)
        @constraint(m, x[1] + x[2] <= 2)
        @constraint(m, x[1] + x[3] <= 2)
        @constraint(m, x[2] + x[3] <= 2)
        @objective(m, Max, x[1] + x[2] - 2*x[3])

        sol = solveMIP(m)

        @test m.ext[:status] == :Optimal
        @test getobjectivevalue(m) == 2
    end
end

# teste PL infeasible
# adicionado por Raphael Saavedra
function test_PL_Infeasible_Raphael(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste PL infeasible" begin
        m = Model(solver = solver)
        @variable(m, x[i=1:2] >= 0)
        @constraint(m, x[1] + x[2] <= -1)
        @objective(m, Max, x[1] + x[2])

        sol = solveMIP(m)

        @test m.ext[:status] == :Infeasible
    end
end

# teste Minimal Unit Commitment
# adicionado por Raphael Saavedra
function test_Minimal_UC(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
  #------------------------------------------------------------------------------
  # Parameters
  T = collect(1:3) # periods
  N = collect(1:3) # generators
  p0 = zeros(N[end]) # initial power output
  v0 = zeros(N[end]) # initial on/off state
  D = [100; 200; 300] # demand
  Cf = [100; 100; 100] # fixed cost
  Cl = [10; 30; 50] # linear cost
  Pmax = [100; 150; 200] # generator capacity
  Pmin = [10; 10; 10] # minimum power output
  RD = [30; 50; 70] # ramp-down limit
  RU = [30; 50; 70] # ramp-up limit
  SD = RD # shutdown ramp limit
  SU = RU # startup ramp limit
  #------------------------------------------------------------------------------
  @testset "Teste Minimal Unit Commitment" begin

    m = Model(solver = solver)
    #----------------------------------------------------------------------------
    # Variables
    @variable(m, p[1:N[end], 1:T[end]] >= 0) # power output
    @variable(m, v[1:N[end], 1:T[end]], Bin) # 1 if generator is on, 0 otherwise
    #----------------------------------------------------------------------------
    # Constraints
    @constraint(m, [t in T], sum(p[i,t] for i in N) == D[t])
    @constraint(m, [t in T, i in N], p[i,t] >= Pmin[i]*v[i,t])
    @constraint(m, [t in T, i in N], p[i,t] <= Pmax[i]*v[i,t])
    @constraint(m, [t in 2:T[end], i in N], p[i,t-1] <= p[i,t] + RD[i]*v[i,t] + SD[i]*(v[i,t-1]-v[i,t]) + Pmax[i]*(1-v[i,t-1]))
    @constraint(m, [t in 2:T[end], i in N], p[i,t] <= p[i,t-1] + RU[i]*v[i,t-1] + SU[i]*(v[i,t]-v[i,t-1]) + Pmax[i]*(1-v[i,t]))
    # Constraints regarding initial state
    @constraint(m, [t in 1, i in N], p0[i] <= p[i,t] + RD[i]*v[i,t] + SD[i]*(v0[i]-v[i,t]) + Pmax[i]*(1-v0[i]))
    @constraint(m, [t in 1, i in N], p[i,t] <= p0[i] + RU[i]*v0[i] + SU[i]*(v[i,t]-v0[i]) + Pmax[i]*(1-v[i,t]))
    #----------------------------------------------------------------------------
    # Objective function
    @objective(m, Min, sum(Cf[i]*v[i,t] for i in N, t in T) + sum(Cl[i]*p[i,t] for i in N, t in T))
  #------------------------------------------------------------------------------
    sol = solveMIP(m)

    @test m.ext[:status] == :Optimal
    @test getobjectivevalue(m) ≈ 17700 atol = 1e-5
    @test getvalue(p) ≈ [30 60 90 ; 50 100 150 ; 20 40 60] atol = 1e-5

  end
end

# teste Sudoku 4x4
# adicionado por Raphael Saavedra
function testSudoku4x4(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste Sudoku 4x4" begin
        n = 4
        m = Model(solver = solver)
        @variable(m, x[i in 1:n, j in 1:n, k in 1:n], Bin)

        fixas = [(1,1,1), (2,2,3), (1,3,4), (3,3,2), (4,4,4), (4,2,1)]
        for idx in fixas
            @constraint(m, x[idx...] == 1)
        end
        @constraint(m, [j in 1:n, k in 1:n], sum(x[:,j,k]) == 1)
        @constraint(m, [i in 1:n, k in 1:n], sum(x[i,:,k]) == 1)
        @constraint(m, [i in 1:n, j in 1:n], sum(x[i,j,:]) == 1)
        @constraint(m, [p in [0,2], q in [0,2], k in 1:n], sum(sum(x[i+p,j+q,k] for i in 1:2) for j in 1:2) == 1)
        @objective(m, Min, 0)

        sol = solveMIP(m)

        M = Matrix(4,4)
        for i = 1 : 4
          for j = 1 : 4
            M[i,j] = find(getvalue(x[i,j,:]).>0)[1]
          end
        end

        @test M == [1 2 4 3; 4 3 1 2; 3 4 2 1; 2 1 3 4]
    end
end


#adicionado por Rodrigo Villas
function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Custo fixo" begin

        #Custo Unitário
        c = [2; 3; 2; 5; 6]


        #Custo Fixo
        f = [1; 3; 1; 5; 10]

        m= Model(solver = solver)

        @variable(m, x[i=1:5]>=0)

        @variable(m,y[i=1:5], Bin)

        @constraint(m, sum(x[j] for j=1:5) >=10)

        @constraint(m,x[1]<=5*y[1])

        @constraint(m,x[2]<=4*y[2])

        @constraint(m,x[3]<=3*y[3])

        @constraint(m,x[4]<=2*y[4])

        @constraint(m,x[5]<=1*y[5])

        @objective(m, Min, sum(f[j]*y[j]+c[j]*x[j] for j=1:5))

        sol = solveMIP(m)
        @test getobjectivevalue(m) == 27
        @test getvalue(x) == [5 2 3 0 0]
end

#-------------------------------------
#adicionado por Rodrigo Villas
function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Cobertura de pontos" begin

        pontosextra= 0

        #um ponto não é coberto por nenhum subconjunto

        S1=[1 0 0 1 ones(1,pontosextra)]
        S2=[1 1 0 0 ones(1,pontosextra)]
        S3=[0 1 0 1 ones(1,pontosextra)]

        c=[4 3 2]

        A=[S1' S2' S3']

        m = Model(solver = solver)
        @variable(m, x[i=1:3], Bin)
        @constraints(m, begin
          constrain[i=1:4+pontosextra], sum(A[i,j]*x[j] for j=1:3)>= 1
          end)
        @objective(m, Min, sum(c[j]*x[j] for j=1:3))

        m.ext[:status] == :Infeasible
end
#-------------------------------------
#adicionado por Rodrigo Villas
function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Cobertura de pontos" begin

        pontosextra= 50

        S1=[1 0 1 1 ones(1,pontosextra)]
        S2=[1 1 0 0 ones(1,pontosextra)]
        S3=[0 1 0 1 ones(1,pontosextra)]

        c=[4 3 2]

        A=[S1' S2' S3']

        m = Model(solver = solver)
        @variable(m, x[i=1:3], Bin)
        @constraints(m, begin
          constrain[i=1:4+pontosextra], sum(A[i,j]*x[j] for j=1:3)>= 1
          end)
        @objective(m, Min, sum(c[j]*x[j] for j=1:3))

        sol = solveMIP(m)
        @test getobjectivevalue(m) == 6
        @test getvalue(x) == [1 0 1]

end


#---------------------

#adicionado por Rodrigo Villas
#Produção com custo fixo Inviavel
function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Custo fixo" begin


        vari=50


        m= Model(solver = solver)

        @variable(m, x[i=1:vari]>=0)

        @variable(m,y[i=1:vari], Bin)

        @constraint(m, sum(x[j] for j=1:vari) >= 2*vari)

        @constraints(m, begin
          constrain[i=1:vari], x[i] <= 1*y[i]
          end)

        @objective(m, Min, sum(j*y[j]+(vari-j)*x[j] for j=1:vari))


        sol = solveMIP(m)
        m.ext[:status] == :Infeasible
end
#------------------------------------------

 #adicionado por Rodrigo Villas
function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Custo fixo" begin

        vari=500

        m= Model(solver = solver)

        @variable(m, x[i=1:vari]>=0)

        @variable(m,y[i=1:vari], Bin)

        @constraint(m, sum(x[j] for j=1:vari) >= vari)

        @constraints(m, begin
          constrain[i=1:vari], x[i] <= 10*y[i]
          end)

        @objective(m, Min, sum(j*y[j]+(vari-j)*x[j] for j=1:vari))

        sol = solveMIP(m)
                @test getobjectivevalue(m) == 36025

                #Os últimos Y's tem que estar "ligados"
end
#-----------------------

#Expansão Unbounded
 #adicionado por Rodrigo Villas
function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Expansão Unbouded" begin


        Cinv = 13.16
           M = 200

               m = Model(solver = solver)
               @variable(m, x[i=1:2])
               @variable(m, u, Bin)
               @objective(m, Min, 4*x[1] + 3*x[2] - u*Cinv)

               @constraint(m, 2*x[1] + 1*x[2] <= 4 +u*M)
               @constraint(m, 1*x[1] + 2*x[2] <= 4 +u*M)

               @constraint(m, 1*x[1] + 0.1*x[2] <= 4 +(1-u)*M)
               @constraint(m, 0.4*x[1] + 1*x[2] <= 4 +(1-u)*M)
        sol = solveMIP(m)
        m.ext[:status] == :Infeasible  #Unbounded
end

#--------------------------
 #adicionado por Rodrigo Villas
function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "PL da minha cabeça" begin


                @variable(m, x[i=1:4]>=0)

                @constraint(m, x[1]+x[2]+x[3]<=3)
                @constraint(m, x[4]+2*x[1]+6*x[3]<=10)
                @constraint(m, 4*x[3]+x[1]+3*x[2]<=5)

                @objective(m, Max, 4*x[1]+5*x[2]+2*x[3]-3*x[4])

                sol = solveMIP(m)
                    @test getobjectivevalue(m) == 13
                    @test getvalue(x) == [2 1 0 0]

end



#--------------------------
 #adicionado por Rodrigo Villas
function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "PL da minha cabeça" begin


                @variable(m, x[i=1:4]>=0)

                @constraint(m, x[1]+x[2]+x[3]<=3)
                @constraint(m, x[4]+2*x[1]+6*x[3]<=10)
                @constraint(m, 4*x[3]+x[1]+3*x[2]<=5)
                @constraint(m, x[3]==2)
                @objective(m, Max, 4*x[1]+5*x[2]+2*x[3]-3*x[4])

                sol = solveMIP(m)
                    m.ext[:status] == :Infeasible

end

function test_optimal_dispatch(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
  @testset "Despacho ótimo" begin
        # Conjuntos

        T = collect(1:2) # Periodos
        L =collect(1:2) # conjunto de linhas
        N = collect(1:2) # conjunto de barras
        I = collect(1:5) # conjunto de usinas
        C = collect(1:7) # conjunto de contingencias

        # Constantes

        x = 1 # reatancia
        c = [10;20;50;200;300] # custo de geração
        u = [10;15;5;20;25] # custo da reserva de subida
        d = [10;15;5;1;1] # custo da reserva de descida
        CU = [100;100;50;10;10] # custo fixo de acoplamento
        CD = CU # custo fixo de desacoplamento
        Pmin = [50;40;10;5;1] # Produção minima de cada usinas
        Pmax = [80;70;60;60;60] # produção maxima de cada usina
        RU = [5;15;25;50;60] # Rampa de subida
        RD = [5;15;50;50;60] # Rampa de descida
        F = [50;50] # capacidade de cada linha
        K = 1 # criterio de segurança

        # Lei de Kirshoff

        A = [1  1;
            -1 -1] # Matriz de incidencia

        S = (A*(1/x)*diagm(ones(length(L))))'

        # Balanço de potencia

        B = [1 1 0 0 1;
             0 0 1 1 0]

        # Demanda

        D = [50;50;50;60;70;80;81;82;83;84;85;80;75;70;80;90;100;100;100;90;80;70;60;50]
        Db=0.5*(D.*ones(24,2))'
        # Contingencias

        a = ones(7,7)-diagm(ones(7))

        # Estado inicial do sistema

        p0 = zeros(5)
        v0 = zeros(5)

        #------------------------------------------------------------------------------
        # Formulação do modelo
        #------------------------------------------------------------------------------

        m = Model(solver = GurobiSolver())

        # Variáveis

        @variable(m, p[1:length(I),1:length(T)] >= 0)
        @variable(m, v[1:length(I),1:length(T)], Bin)
        @variable(m, f[1:length(L),1:length(T)])
        @variable(m, θ[1:length(N),1:length(T)])
        @variable(m, ru[1:length(I),1:length(T)] >= 0)
        @variable(m, rd[1:length(I),1:length(T)] >= 0)
        @variable(m, pk[1:length(I),1:length(T),1:7] >= 0)
        @variable(m, fk[1:length(L),1:length(T),1:7])
        @variable(m, θk[1:length(N),1:length(T),1:7])
        @variable(m, 0 <= vu[1:length(I),1:length(T)] <= 1)
        @variable(m, 0 <= vd[1:length(I),1:length(T)] <= 1)

        # Função objetivo

        @objective(m, Min,sum(c[i]*p[i,t]+CU[i]*vu[i,t]+CD[i]*vd[i,t]+u[i]*ru[i,t]+d[i]*rd[i,t] for i in I, t in T))

        # Restições

        @constraint(m, [i in I,t in 2:T[end]], vu[i,t] - vd[i,t] == v[i,t] - v[i,t-1])
        @constraint(m, [i in I,t in T], vu[i,t] <= v[i,t])
        @constraint(m, [i in I,t in T], vd[i,t] <= 1- v[i,t])
        @constraint(m, [i in I,t in T], Pmin[i]*v[i,t] <= p[i,t])
        @constraint(m, [i in I,t in T], p[i,t] <= Pmax[i]*v[i,t])
        @constraint(m, [i in I, t in T[2:end]], p[i,t] - p[i,t-1] <= RU[i]*v[i,t-1] + Pmax[i]*vu[i,t])
        @constraint(m, [i in I, t in T[2:end]], p[i,t-1] - p[i,t] <= RD[i]*v[i,t] + Pmax[i]*vd[i,t])
        @constraint(m, pre[b in N,t in T], sum(A[b,l]*f[l,t] for l in L) + sum(B[b,i]*p[i,t] for i in I) == Db[b,t])
        @constraint(m, [l in L,t in T], f[l,t] == sum(S[l,b]*θ[b,t] for b in N))
        @constraint(m, [l in L,t in T], -F[l] <= f[l,t])
        @constraint(m, [l in L,t in T], f[l,t] <= F[l])
        @constraint(m, [i in I,t in T], p[i,t] + ru[i,t] <= Pmax[i]*v[i,t])
        @constraint(m, [i in I,t in T], p[i,t] - rd[i,t] >= Pmin[i]*v[i,t])
        @constraint(m, [i in I,t in T], ru[i,t] <= 0.9*RU[i])
        @constraint(m, [i in I,t in T], rd[i,t] <= 0.9*RD[i])


        if K!=0
          @constraint(m, pos[b in N,t in T,k in C], sum(A[b,l]*fk[l,t,k] for l in L) + sum(B[b,i]*pk[i,t,k] for i in I) == Db[b,t])
          @constraint(m, [l in L,t in T,k in C], fk[l,t,k] == a[k,l+5]*sum(S[l,b]*θk[b,t,k] for b in N))
          @constraint(m, [l in L,t in T,k in C], -F[l]*a[k,l+5] <= fk[l,t,k])
          @constraint(m, [l in L,t in T,k in C], fk[l,t,k] <= F[l]*a[k,l+5])
          @constraint(m, [i in I,t in T,k in C], a[k,i]*(p[i,t]-rd[i,t]) <= pk[i,t,k])
          @constraint(m, [i in I,t in T,k in C], pk[i,t,k] <= a[k,i]*(p[i,t]+ru[i,t]))
        end


        @constraint(m, [i in I,t in 1], vu[i,t] - vd[i,t] == v[i,t] - v0[i])
        @constraint(m, [i in I, t in 1], p[i,t] - p0[i] <= RU[i]*v0[i] + Pmax[i]*vu[i,t])
        @constraint(m, [i in I, t in 1], p0[i] - p[i,t] <= RD[i]*v[i,t] + Pmax[i]*vd[i,t])
        sol = solveMIP(m)

        @test getobjectivevalue(m) == 5150.0
        @test sum(getvalue(p[2,:])) == 98.0
        @test sum(getvalue(v[1,:])) == 0.00
        @test sum(getvalue(v)) == 4.00
  end
end
