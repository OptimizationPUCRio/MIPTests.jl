# using JuMP, using Base.test

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
