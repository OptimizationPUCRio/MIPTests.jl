include("miptests.jl")

testset = [
    (test1, "JG", "Example"),

    (test2, "GB", "MIP-p"),
    (test_P1_Guilherme, "GB", "MIP-g1"),
    (test_TSPbin7_Guilherme, "GB", "MIP-g1-TSP-bin"),
    (test_TSPmip7_Guilherme, "GB", "MIP-g1-TSP-mip"),
    (test_PL_Infeasible_Guilherme, "GB", "MIP-inf-p"),
    (test_MIP_medio_Guilherme, "GB", "MIP-m"),
    (test_MIP_Grande_Guilherme, "GB", "MIP-g"),
    (test_PL_Unbounded_Guilherme, "GB", "LP-unb"),
    (test_MIP_Unbounded_Guilherme, "GB", "MIP-unb"),
    (test_MIP_Infeasible_Minimal_Guilherme, "GB", "MIP-inf-m"),
    (test_MIP_Infeasible_Pequeno_Guilherme, "GB", "MIP-inf-m"),
    (test_PL_Guilherme, "GB", "LP-p"),

    (test_P1_Brito, "EB", "MIP-g1"),
    (test_PL_Simples_Brito, "EB", "LP-p"),
    (test_PL_Infeasible_Brito, "EB", "LP-inf-p"),
    (test_PL_Unbounded_Brito, "EB", "LP-unb"),
    (test_MIP_Minimal_Brito, "EB", "MIP-pp"),
    (test_MIP_Pequeno_Brito, "EB", "MIP-p"),

    (testSudoku, "RS", "MIP-m"),
    (testInfeasibleKnapsack, "RS", "MIP-inf-p"),
    (testRobustCCUC, "RS", "MIP-g1"),
    (testUnboundedKnapsack, "RS", "MIP-unb"),
    (testInfeasibleUC, "RS", "MIP-inf"),
    (test_PL_Simples_Raphael, "RS", "LP-opt"),
    (test_PL_Infeasible_Raphael, "RS", "LP-inf"),
    (test_Minimal_UC, "RS", "MIP-pp"),
    (testSudoku4x4, "RS", "MPI-p"),

    (testCaminho, "CG", "MIP-p"),
    (test_optimal_dispatch, "CG", "MIP-g"),

    (test3, "AR", "MIP-p"),
    (test3_2, "AR", "MIP-unb"),
    (test3_3, "AR", "MPI-inf"),
    (test_feature_selection_pequeno_viavel, "AR", "MIP-p"),
    (test_feature_selection_medio, "AR", "MIP-m"),
    (test_feature_selection_grande, "AR", "MIP-g"),
    (test_feature_selection_pequeno_inviavel, "AR", "MIP-inf-m"),
    (teste_PL_andrew_unbounded, "AR", "LP-unb"),
    (teste_PL_andrew_viavel, "AR", "LP-opt"),
    (teste_PL_andrew_inviavel, "AR", "LP-inf"),

    (test_P1_Andrew_Bianca_viavel, "AR", "MIP-g1"),

    (test_rv_1, "RV", "MIP-pp"),
    (test_rv_3, "RV", "MIP-p"),
    (test_rv_5, "RV", "MIP-m"),
    (test_rv_2, "RV", "MIP-inf-p"),
    (test_rv_4, "RV", "LP-inf"),
    (test_rv_7, "RV", "LP-opt"),
    (test_rv_8, "RV", "LP-inf"),
    (test_rv_p1, "RV", "MIP-p1"),

]

function runtests(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver();
    author = String[], kind = String[], ignore = String[], name = "")

    table = MIPSolution[]
    @testset "Main" begin
        for teste in testset
            if !("$(teste[1])" in ignore) && (isempty(author) || teste[2] in author) && (isempty(kind) || teste[3] in kind)
                println("$(teste[1])")
                line = teste[1](solveMIP, solver)
                line.name = "$(teste[1]) - $(teste[2])"
                push!(table, line)
            end
        end
        mat = table2mat(table)
        if name == ""
            name = "result_$(solveMIP)_$(solver)"
        end
        writecsv(name*".out",mat)
    end

    return table
end

function table2mat(table)
    out = Any["Name" "Pass" "Objective" "BestBound" "Time" "Nodes" "IntSols" "Status"]
    for line in table
        out = vcat(out, Any[line.name line.pass line.objective line.bestbound line.time line.nodes line.intsols line.status])
    end
    return out
end

function jumpsolve(m)
    tic()
    sol = solve(m)
    time = toq()
    m.ext[:status] = sol
    m.ext[:time] = time
end
