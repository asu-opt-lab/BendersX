using BendersDecomposition
using Test
using JuMP

@testset verbose = true "UFLP Auto Decompose - Sequential" begin
    instances = setdiff(1:71, [67])

    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data
            problem = read_uflp_benchmark_data("p$(i)")


            # Algorithm parameters
            benders_param = BendersSeqParam(
                time_limit = 200.0,
                gap_tolerance = 1e-6,
                verbose = false
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
                @info "solving UFLP p$i - auto_decompose - seq - classical"

                data, master, oracle = auto_decompose(
                    model;
                    master_solver_param = master_solver_param,
                    oracle_solver_param = oracle_solver_param
                )

                env = BendersSeq(data, master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end

            @testset "Disjunctive oracle" begin
                @info "solving UFLP p$i - auto_decompose - seq - disjunctive"

                dcglp_param = DcglpParam(
                    time_limit = 1000.0,
                    gap_tolerance = 1e-3,
                    halt_limit = 3,
                    iter_limit = 250,
                    verbose = false
                )

                oracle_param = DisjunctiveOracleParam(
                    norm = LpNorm(1.0),
                    split_index_selection_rule = RandomFractional(),
                    disjunctive_cut_append_rule = AllDisjunctiveCuts(),
                    strengthened = true,
                    add_benders_cuts_to_master = true,
                    fraction_of_benders_cuts_to_master = 1.0,
                    reuse_dcglp = true,
                    adjust_t_to_fx = true
                )

                dcglp_solver_param = Dict("solver" => "Gurobi", "InfUnbdInfo" => 1)

                data, master, oracle = auto_decompose(
                    model,
                    :disjunctive;
                    master_solver_param = master_solver_param,
                    oracle_param = oracle_param,
                    typical_oracle_solver_param = oracle_solver_param,
                    dcglp_solver_param = dcglp_solver_param,
                    dcglp_param = dcglp_param
                )

                env = BendersSeq(data, master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end
        end
    end
end
