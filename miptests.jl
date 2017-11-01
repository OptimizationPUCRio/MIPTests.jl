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

# teste Sudoku
# adicionado por Raphael Saavedra
function testSudoku(solveMIP::Function)
  @testset "Teste Sudoku" begin
    n = 9
    model = Model()
    @variable(model, x[i in 1:n, j in 1:n, k in 1:n], Bin)
    # ============================================
    # Definir casas preenchidas
    @constraint(model, x[1,3,4] == 1)
    @constraint(model, x[1,5,6] == 1)
    @constraint(model, x[1,9,2] == 1)
    @constraint(model, x[2,1,8] == 1)
    @constraint(model, x[2,3,5] == 1)
    @constraint(model, x[2,6,2] == 1)
    @constraint(model, x[2,8,3] == 1)
    @constraint(model, x[3,5,3] == 1)
    @constraint(model, x[3,8,6] == 1)
    @constraint(model, x[4,2,2] == 1)
    @constraint(model, x[4,3,8] == 1)
    @constraint(model, x[5,6,4] == 1)
    @constraint(model, x[6,1,7] == 1)
    @constraint(model, x[6,5,5] == 1)
    @constraint(model, x[6,9,9] == 1)
    @constraint(model, x[7,3,2] == 1)
    @constraint(model, x[7,6,1] == 1)
    @constraint(model, x[8,2,7] == 1)
    @constraint(model, x[8,5,4] == 1)
    @constraint(model, x[8,7,9] == 1)
    @constraint(model, x[8,8,5] == 1)
    @constraint(model, x[9,1,6] == 1)
    @constraint(model, x[9,8,4] == 1)
    # ============================================
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
