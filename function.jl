

module SolverBrendon

  using JuMP, Gurobi

  mutable struct prob_lp
      A
      c
      n
      m
      xlb
      xub
      l
      u
      vtypes
      solver
  end

  mutable struct resposta_relaxado
      obj
      vars
      status
  end

  mutable struct modelo_lista
      problem::prob_lp
      resp::resposta_relaxado
      modelo_lista() = new()
  end



  function converte_modelo(mw::Model)

    m = deepcopy(mw)

    if m.objSense == :Max
        @objective(m,:Min,-m.obj)
    end
    A = JuMP.prepConstrMatrix(m)
    n, mm= size(A)
    xlb = copy(m.colLower)
    xub = copy(m.colUpper)
    l,u = JuMP.prepConstrBounds(m)
    c = JuMP.prepAffObjective(m)
    solver = m.solver
    vtypes = copy(m.colCat)
    problem = prob_lp(A,c,n,mm,xlb,xub,l,u,vtypes,solver)


    return problem
  end


  function solve_relax(problema::prob_lp)


    mod = Model(solver=problema.solver)
    @variable(mod, x[1:problema.m])
    for i in 1:problema.m
        setlowerbound(x[i], problema.xlb[i])
        setupperbound(x[i], problema.xub[i])

    end
    @constraint(mod, problema.l .<= problema.A*x .<= problema.u)
    @objective(mod, Min, dot(problema.c,x))

    status = solve(mod)
    resp = resposta_relaxado(getobjectivevalue(mod),getvalue(x),status)
    return resp
  end

  function podas(resp ::resposta_relaxado, global_bound, v_bin )

    #Poda por Inviabilidade###########################
    if(resp.status != :Optimal)
      return "erro"
    end
    #Poda por Limite#####################################################
    if (resp.obj>global_bound[2])
      return "erro"
    end
    #Poda por Otimalidade###############################################
    if ( sum(abs.(resp.vars - round.(resp.vars)).*v_bin) == 0 )
      return Float64(resp.obj)
    end
    return "sucesso"
  end

  function exp_mod(m,zinf,nos_explorados,int_sol,time,global_bound)

    m.objVal = zinf.resp.obj
    m.colVal = zinf.resp.vars
    m.objBound =  global_bound


    m.ext[:nos] = nos_explorados
    m.ext[:solucao_inteira] = int_sol
    m.ext[:status] = zinf.resp.status
    m.ext[:time] = time

    return m
  end

  function SolveMIP(m)

    # Inicia o Clock###############################################
    tic()
    # Criacao da lista: ###########################################
    lista = Vector{modelo_lista}()
    nos_explorados = Vector{modelo_lista}()
    int_sol = Vector{modelo_lista}()
    # Adicionando primeiro problema, de forma manual: #############
    zinf = modelo_lista()
    zinf.problem = converte_modelo(m)
    zinf.resp = solve_relax(zinf.problem)
    global_bound = [zinf.resp.obj,+Inf]
    push!(lista,zinf)

    #vetor vtype###################################################
    v_bin=zeros(size(zinf.problem.vtypes)[1])

    for i = 1: size(v_bin)[1]
      if((zinf.problem.vtypes[i]==:Int)|(zinf.problem.vtypes[i]==:Bin))
        v_bin[i]=1
      end
    end


    if(sum(v_bin)!=0)
      ############################### ###############################
      iter = 0
      while (abs(global_bound[2] - global_bound[1]) >= 0.00000005 && size(lista)[1] != 0)
          #Seleciona problema  ##########################################
          # ind_prob = #ind do problema selecionado
          ind_prob = 1

          ############################### ###############################

          #Select variables #############################################


          ind = indmax(abs.(lista[ind_prob].resp.vars.*v_bin - round.(lista[ind_prob].resp.vars.*v_bin)))
          ############################### ###############################

          #Branch #######################################################
          prob_lb = deepcopy(lista[ind_prob].problem)
          prob_lb.xlb[ind] = floor(lista[ind_prob].resp.vars[ind]) + 1

          prob_ub = deepcopy(lista[ind_prob].problem)
          prob_ub.xub[ind] = floor(lista[ind_prob].resp.vars[ind])
          ############################### ###############################

          #Solve das folhas #############################################
          resp_lb = solve_relax(prob_lb)
          resp_ub = solve_relax(prob_ub)
          ############################### ###############################

          #Podas ########################################################
          poda_lb = podas(resp_lb, global_bound,v_bin)
          poda_ub = podas(resp_ub, global_bound,v_bin)
          ############################### ###############################

          #Monta problema  ##############################################
          novo_lb = 0
          if typeof(poda_lb) == Float64
              global_bound[2] = poda_lb
              novo_lb = modelo_lista()
              novo_lb.problem = prob_lb
              novo_lb.resp = resp_lb
              push!(int_sol,novo_lb)
          elseif poda_lb == "sucesso"
              novo_lb = modelo_lista()
              novo_lb.problem = prob_lb
              novo_lb.resp = resp_lb
          end

          novo_ub = 0
          if (typeof(poda_ub) == Float64)
              if (typeof(poda_lb) == Float64)
                  if (poda_ub < poda_lb)
                      global_bound[2] = poda_ub
                      novo_ub = modelo_lista()
                      novo_ub.problem = prob_ub
                      novo_ub.resp = resp_ub
                      push!(int_sol,novo_ub)
                  else
                    global_bound[2] = poda_ub
                    novo_ub = modelo_lista()
                    novo_ub.problem = prob_ub
                    novo_ub.resp = resp_ub
                    push!(int_sol,novo_ub)
                end

              else
                  global_bound[2] = poda_ub
                  novo_ub = modelo_lista()
                  novo_ub.problem = prob_ub
                  novo_ub.resp = resp_ub
                  push!(int_sol,novo_ub)
              end
          elseif poda_ub == "sucesso"
              novo_ub = modelo_lista()
              novo_ub.problem = prob_ub
              novo_ub.resp = resp_ub
          end
          ############################### ###############################

          #Atualiza Zinf ################################################
          if typeof(poda_lb) == Float64
              if novo_lb.resp.obj <= global_bound[2]
                  zinf = deepcopy(novo_lb)
              end
          end
          if typeof(poda_ub) == Float64
              if typeof(poda_lb) == Float64
                  if novo_ub.resp.obj <= global_bound[2] && novo_ub.resp.obj <= novo_lb.resp.obj
                      zinf = deepcopy(novo_ub)
                  end
              else
                  if novo_ub.resp.obj <= global_bound[2]
                      zinf = deepcopy(novo_ub)
                  end
              end
          end
          ############################### ###############################

          #Remove o problema original e adiciona os novos ###############
          #Remove

          push!(nos_explorados, lista[ind_prob])
          deleteat!(lista,ind_prob)

          #adiciona
          if poda_lb == "sucesso"
              push!(lista,novo_lb)
          end
          if poda_ub == "sucesso"
              push!(lista,novo_ub)
          end
          iter += 1
      end
    end
    if m.objSense == :Max
        zinf.resp.obj = - zinf.resp.obj
    end

    time = toc()

    model = exp_mod(m,zinf,nos_explorados,int_sol,time,global_bound)

    return model
  end
end
