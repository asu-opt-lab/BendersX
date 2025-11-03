using BendersDecomposition
using Test
using JuMP

@testset verbose = true "CFLP Auto Decompose - Sequential" begin
    instances = setdiff(1:71, [67])

    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data
            problem = read_cflp_benchmark_data("p$i")

            # Create traditional data for MIP reference
            dim_x = problem.n_facilities
            dim_t = 1
            c_x = problem.fixed_costs
            c_t = [1.0]
            data_for_mip = Data(dim_x, dim_t, problem, c_x, c_t)

            # Algorithm parameters
            benders_param = BendersSeqParam(
                time_limit = 200.0,
                gap_tolerance = 1e-6,
                verbose = false
            )

            # Solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPXPARAM_Threads" => 4)
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9, "CPXPARAM_Threads" => 4)
            oracle_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)

            # Solve MIP for reference
            mip = Mip(data_for_mip)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data_for_mip)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)

            @testset "Classical oracle" begin
                @info "solving CFLP p$i - auto_decompose - seq - classical"

                data, master, oracle = auto_decompose(
                    mip.model;
                    master_solver_param = master_solver_param,
                    oracle_solver_param = oracle_solver_param
                )

                env = BendersSeq(data_for_mip, master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end

            @testset "Disjunctive oracle" begin
                @info "solving CFLP p$i - auto_decompose - seq - disjunctive"

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

                dcglp_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)

                data, master, oracle = auto_decompose(
                    mip.model,
                    :disjunctive;
                    master_solver_param = master_solver_param,
                    oracle_param = oracle_param,
                    typical_oracle_solver_param = oracle_solver_param,
                    dcglp_solver_param = dcglp_solver_param,
                    dcglp_param = dcglp_param
                )

                env = BendersSeq(data_for_mip, master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end
        end
    end
end
