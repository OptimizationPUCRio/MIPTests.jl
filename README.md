# MIPTests.jl

A collection of MIP instances

## Usage

The main script is **runmiptests.jl**

In order to test a MIP solver in this testbench use the function:

`runtests(solvemethod, solver)`

Where `solver` is a typical JuMP solver like: `GurobiSolver(OutputFlag = 0)`, `Xpress.XpressSolver()`, `CplexSolver(CPX_PARAM_SCRIND = 0)`...

And `solvemethod` is a julia solve functions. To use JuMP's solver default use: `jumpsolve` defined in this script

## Example:

```
include("runmiptests.jl")
using Gurobi
runtests(jumpsolve, GurobiSolver(OutputFlag = 0))
```

This will write the file: `result_jumpsolve_Grobi.GurobiSolver(nothing, Any[(:OutputFlag ,0)]).out`
