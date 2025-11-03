using BendersDecomposition
using Test
using JuMP

@testset verbose = true "UFLP Auto Decompose - Callback" begin
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

            # Solve MIP for reference
            mip = Mip(data_for_mip)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data_for_mip)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)

            @testset "Classical oracle" begin
                @info "solving UFLP p$i - auto_decompose - callback - classical"

                data, master, oracle = auto_decompose(
                    mip.model;
                    master_solver_param = master_solver_param,
                    oracle_solver_param = oracle_solver_param
                )

                root_seq_type = BendersSeq
                root_param = BendersSeqParam(
                    time_limit = 200.0,
                    gap_tolerance = 1e-6,
                    verbose = false
                )
                root_preprocessing = RootNodePreprocessing(oracle, root_seq_type, root_param)
                lazy_callback = LazyCallback(oracle)
                user_callback = NoUserCallback()

                env = BendersBnB(data_for_mip, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end

            @testset "Disjunctive oracle" begin
                # First auto_decompose: classical oracle for lazy callback
                data_classical, master_classical, lazy_oracle_classical = auto_decompose(
                    mip.model;
                    master_solver_param = master_solver_param,
                    oracle_solver_param = oracle_solver_param
                )

                # Second auto_decompose: disjunctive oracle for user callback
                data_disj, master_disj, disjunctive_oracle = auto_decompose(
                    mip.model,
                    :disjunctive;
                    master_solver_param = master_solver_param,
                    typical_oracle_solver_param = oracle_solver_param,
                    dcglp_solver_param = dcglp_solver_param,
                    dcglp_param = dcglp_param
                )

                # Test various parameter combinations
                for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [true], lift in [true], p in [1.0], disjunctive_cut_append_rule in [AllDisjunctiveCuts()], adjust_t_to_fx in [true]

                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; lift $lift; p $p; dcut_append $disjunctive_cut_append_rule; adjust_t_to_fx $adjust_t_to_fx" begin
                        # Configure DisjunctiveOracleParam
                        oracle_param = DisjunctiveOracleParam(
                            norm = LpNorm(p),
                            split_index_selection_rule = RandomFractional(),
                            disjunctive_cut_append_rule = disjunctive_cut_append_rule,
                            strengthened = strengthened,
                            add_benders_cuts_to_master = add_benders_cuts_to_master,
                            fraction_of_benders_cuts_to_master = 0.5,
                            reuse_dcglp = reuse_dcglp,
                            lift = lift,
                            adjust_t_to_fx = adjust_t_to_fx
                        )
                        set_parameter!(disjunctive_oracle, oracle_param)
                        update_model!(disjunctive_oracle, data_for_mip)

                        @testset "NoSeq" begin
                            @info "solving UFLP p$i - auto_decompose - callback - disjunctive - no seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule adjust_t_to_fx $adjust_t_to_fx"
                            master = Master(data_for_mip; solver_param = master_solver_param)
                            update_model!(master, data_for_mip)
                            root_preprocessing = NoRootNodePreprocessing()
                            lazy_callback = LazyCallback(lazy_oracle_classical)
                            user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
                            env = BendersBnB(data_for_mip, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                        end

                        @testset "Seq" begin
                            @info "solving UFLP p$i - auto_decompose - callback - disjunctive - seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule adjust_t_to_fx $adjust_t_to_fx"
                            master = Master(data_for_mip; solver_param = master_solver_param)
                            update_model!(master, data_for_mip)
                            root_preprocessing = RootNodePreprocessing(lazy_oracle_classical, BendersSeq, BendersSeqParam(time_limit=200.0, gap_tolerance=1e-6, verbose=false))
                            lazy_callback = LazyCallback(lazy_oracle_classical)
                            user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
                            env = BendersBnB(data_for_mip, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                        end

                        @testset "SeqInOut" begin
                            @info "solving UFLP p$i - auto_decompose - callback - disjunctive - seqinout - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule adjust_t_to_fx $adjust_t_to_fx"
                            master = Master(data_for_mip; solver_param = master_solver_param)
                            update_model!(master, data_for_mip)
                            root_preprocessing = RootNodePreprocessing(lazy_oracle_classical, BendersSeqInOut, BendersSeqInOutParam(time_limit = 300.0, gap_tolerance = 1e-6, stabilizing_x = ones(data_for_mip.dim_x), α = 0.9, λ = 0.1, verbose = false))
                            lazy_callback = LazyCallback(lazy_oracle_classical)
                            user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
                            env = BendersBnB(data_for_mip, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                        end
                    end
                end
            end
        end
    end
end
