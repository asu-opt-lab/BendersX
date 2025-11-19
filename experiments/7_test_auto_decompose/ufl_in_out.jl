using BendersDecomposition
using Test
using JuMP

@testset verbose = true "UFLP Auto Decompose - Sequential In/Out" begin
    instances = setdiff(1:71, [67])

    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data
            problem = read_uflp_benchmark_data("p$(i)")

            # Algorithm parameters
            benders_inout_param = BendersSeqInOutParam(
                time_limit = 200.0,
                gap_tolerance = 1e-6,
                verbose = false,
                stabilizing_x = ones(problem.n_facilities),
                α = 0.9,
                λ = 0.1
            )

            # Solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPXPARAM_Threads" => 4)
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9, "CPXPARAM_Threads" => 4)
            oracle_solver_param = Dict("solver" => "Gurobi", "InfUnbdInfo" => 1)

            # Solve MIP for reference
            model = Model()
            set_optimizer_attribute(model, MOI.Silent(), true)
            assign_attributes!(model, mip_solver_param)
            N, M = problem.n_facilities, problem.n_customers
            @variable(model, x[1:N], Bin)
            @variable(model, y[1:N, 1:M] >= 0)
            @objective(model, Min, 
                sum(problem.costs[i,j] * problem.demands[j] * y[i,j] for i in 1:N, j in 1:M) + 
                sum(problem.fixed_costs[i] * x[i] for i in 1:N)
            )
            @constraint(model, demand[j in 1:M], sum(y[:,j]) == 1)
            @constraint(model, facility_open[i in 1:N, j in 1:M], y[i,j] <= x[i])
            optimize!(model)
            @assert termination_status(model) == OPTIMAL
            mip_opt_val = objective_value(model)

            @testset "Classical oracle" begin
                @info "solving UFLP p$i - auto_decompose - seqInOut - classical"

                data, master, oracle = auto_decompose(
                    model;
                    master_solver_param = master_solver_param,
                    oracle_solver_param = oracle_solver_param
                )

                env = BendersSeqInOut(data, master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end
        end
    end
end
