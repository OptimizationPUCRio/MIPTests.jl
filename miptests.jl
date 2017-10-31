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
