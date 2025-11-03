using BendersDecomposition
using Test
using JuMP

@testset verbose = true "UFLP Auto Decompose - Sequential In/Out" begin
    instances = setdiff(1:71, [67])

    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data
            problem = read_uflp_benchmark_data("p$(i)")

            # Create traditional data for MIP reference
            dim_x = problem.n_facilities
            dim_t = 1
            c_x = problem.fixed_costs
            c_t = [1.0]
            data_for_mip = Data(dim_x, dim_t, problem, c_x, c_t)

            # Algorithm parameters
            benders_inout_param = BendersSeqInOutParam(
                time_limit = 200.0,
                gap_tolerance = 1e-6,
                verbose = false,
                stabilizing_x = ones(data_for_mip.dim_x),
                α = 0.9,
                λ = 0.1
            )

            # Solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPXPARAM_Threads" => 4)
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9, "CPXPARAM_Threads" => 4)
            oracle_solver_param = Dict("solver" => "Gurobi", "InfUnbdInfo" => 1)

            # Solve MIP for reference
            mip = Mip(data_for_mip)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data_for_mip)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)

            @testset "Classical oracle" begin
                @info "solving UFLP p$i - auto_decompose - seqInOut - classical"

                data, master, oracle = auto_decompose(
                    mip.model;
                    master_solver_param = master_solver_param,
                    oracle_solver_param = oracle_solver_param
                )

                env = BendersSeqInOut(data_for_mip, master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end
        end
    end
end
