using BendersDecomposition
using Test
using JuMP

# Algorithm parameters
benders_param = BendersBnBParam(
    time_limit = 200.0,
    gap_tolerance = 1e-6,
    verbose = false
)
# Solver parameters
mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPXPARAM_Threads" => 4)
master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9, "CPXPARAM_Threads" => 4)
oracle_solver_param = Dict("solver" => "Gurobi", "InfUnbdInfo" => 1, "FeasibilityTol" => 1e-9, "OptimalityTol" => 1e-9)
dcglp_solver_param = Dict("solver" => "Gurobi", "InfUnbdInfo" => 1, "FeasibilityTol" => 1e-9, "OptimalityTol" => 1e-9)
dcglp_param = DcglpParam(
    time_limit = 200.0,
    gap_tolerance = 1e-3,
    halt_limit = 3,
    iter_limit = 15,
    verbose = false
)
user_cb_param = UserCallbackParam(frequency=10)

@testset verbose = true "UFLP Auto Decompose - Callback" begin
    instances = setdiff(1:71, [67])

    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data
            problem = read_uflp_benchmark_data("p$i")

            # Create traditional data for MIP reference
            dim_x = problem.n_facilities
            dim_t = 1
            c_x = problem.fixed_costs
            c_t = [1.0]
            data = Data(dim_x, dim_t, problem, c_x, c_t)

            # Solve MIP for reference
            mip = Mip(data)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)

            @testset "Classical oracle" begin
                @testset "NoSeq" begin
                    @info "solving UFLP p$i - classical oracle - no seq..."
                    data, master, typical_oracle = auto_decompose(
                        mip.model;
                        master_solver_param = master_solver_param,
                        oracle_solver_param = oracle_solver_param
                    )
                    root_preprocessing = NoRootNodePreprocessing()
                    lazy_callback = LazyCallback(typical_oracle)
                    user_callback = NoUserCallback()
                    env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end
                @testset "Seq" begin
                    @info "solving UFLP p$i - classical oracle - seq..."
                    data, master, typical_oracle = auto_decompose(
                        mip.model;
                        master_solver_param = master_solver_param,
                        oracle_solver_param = oracle_solver_param
                    )
                    root_seq_type = BendersSeq
                    root_param = BendersSeqParam(;
                        time_limit = 200.0,
                        gap_tolerance = 1e-6,
                        verbose = false
                    )
                    root_preprocessing = RootNodePreprocessing(typical_oracle, root_seq_type, root_param)
                    lazy_callback = LazyCallback(typical_oracle)
                    user_callback = NoUserCallback()
                    env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end
                @testset "SeqInOut" begin
                    @info "solving UFLP p$i - classical oracle - seqinout..."
                    data, master, typical_oracle = auto_decompose(
                        mip.model;
                        master_solver_param = master_solver_param,
                        oracle_solver_param = oracle_solver_param
                    )
                    root_seq_type = BendersSeqInOut
                    root_param = BendersSeqInOutParam(
                        time_limit = 300.0,
                        gap_tolerance = 1e-6,
                        stabilizing_x = ones(data.dim_x),
                        α = 0.9,
                        λ = 0.1,
                        verbose = false
                    )
                    root_preprocessing = RootNodePreprocessing(typical_oracle, root_seq_type, root_param)
                    lazy_callback = LazyCallback(typical_oracle)
                    user_callback = NoUserCallback()
                    env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end
            end

            @testset "Disjunctive oracle" begin

                oracle_param = DisjunctiveOracleParam(
                    norm = LpNorm(1.0),
                    split_index_selection_rule = RandomFractional(),
                    disjunctive_cut_append_rule = AllDisjunctiveCuts(),
                    strengthened = true,
                    add_benders_cuts_to_master = true,
                    fraction_of_benders_cuts_to_master = 1.0,
                    reuse_dcglp = true,
                    adjust_t_to_fx = false
                )

                @testset "NoSeq" begin
                    @info "solving UFLP p$i - disjunctive oracle - no seq..."
                    data, master, disjunctive_oracle = auto_decompose(
                        mip.model,
                        :disjunctive;
                        master_solver_param = master_solver_param,
                        oracle_param = oracle_param,
                        typical_oracle_solver_param = oracle_solver_param,
                        dcglp_solver_param = dcglp_solver_param,
                        dcglp_param = dcglp_param
                    )
                    _, _, lazy_oracle = auto_decompose(
                        mip.model;
                        master_solver_param = master_solver_param,
                        oracle_solver_param = oracle_solver_param
                    )
                    root_preprocessing = NoRootNodePreprocessing()
                    lazy_callback = LazyCallback(lazy_oracle)
                    user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
                    env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end
                @testset "Seq" begin
                    @info "solving UFLP p$i - disjunctive oracle - seq..."
                    data, master, disjunctive_oracle = auto_decompose(
                        mip.model,
                        :disjunctive;
                        master_solver_param = master_solver_param,
                        oracle_param = oracle_param,
                        typical_oracle_solver_param = oracle_solver_param,
                        dcglp_solver_param = dcglp_solver_param,
                        dcglp_param = dcglp_param
                    )
                    _, _, lazy_oracle = auto_decompose(
                        mip.model;
                        master_solver_param = master_solver_param,
                        oracle_solver_param = oracle_solver_param
                    )
                    root_seq_type = BendersSeq
                    root_param = BendersSeqParam(;
                        time_limit = 200.0,
                        gap_tolerance = 1e-6,
                        verbose = false
                    )
                    root_preprocessing = RootNodePreprocessing(lazy_oracle, root_seq_type, root_param)
                    lazy_callback = LazyCallback(lazy_oracle)
                    user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
                    env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end
                @testset "SeqInOut" begin
                    @info "solving UFLP p$i - disjunctive oracle - seqinout..."
                    data, master, disjunctive_oracle = auto_decompose(
                        mip.model,
                        :disjunctive;
                        master_solver_param = master_solver_param,
                        oracle_param = oracle_param,
                        typical_oracle_solver_param = oracle_solver_param,
                        dcglp_solver_param = dcglp_solver_param,
                        dcglp_param = dcglp_param
                    )
                    _, _, lazy_oracle = auto_decompose(
                        mip.model;
                        master_solver_param = master_solver_param,
                        oracle_solver_param = oracle_solver_param
                    )
                    root_seq_type = BendersSeqInOut
                    root_param = BendersSeqInOutParam(
                        time_limit = 300.0,
                        gap_tolerance = 1e-6,
                        stabilizing_x = ones(data.dim_x),
                        α = 0.9,
                        λ = 0.1,
                        verbose = false
                    )
                    root_preprocessing = RootNodePreprocessing(lazy_oracle, root_seq_type, root_param)
                    lazy_callback = LazyCallback(lazy_oracle)
                    user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
                    env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end
            end
        end
    end
end
