using Test
using BendersLibrary
using BendersBase
using JuMP

@testset "Auto Decompose Disjunctive" begin

    @testset "Basic Functionality" begin
        @testset "DisjunctiveOracle basic structure" begin
            # Build a simple facility location model
            model = Model()
            @variable(model, x[1:3], Bin)
            @variable(model, y[1:3, 1:2] >= 0)
            @objective(model, Min, 10*x[1] + 15*x[2] + 12*x[3] + sum(y))
            @constraint(model, sum(x) >= 1)
            @constraint(model, demand[j in 1:2], sum(y[:,j]) == 1)
            @constraint(model, capacity[i in 1:3, j in 1:2], y[i,j] <= x[i])

            # Decompose the model with disjunctive oracle
            data, master, oracle = auto_decompose(model, :disjunctive)

            # Verify Data structure
            @test data.dim_x == 3
            @test data.dim_t == 1
            @test length(data.c_x) == 3
            @test length(data.c_t) == 1

            # Verify Master model structure
            @test num_variables(master.model) == 4  # 3 x variables + 1 t variable
            @test length(master.x) == 3
            @test length(master.t) == 1

            # Verify DisjunctiveOracle structure
            @test oracle isa DisjunctiveOracle
            @test length(oracle.typical_oracles) == 2
            @test oracle.typical_oracles[1] isa ClassicalOracle
            @test oracle.typical_oracles[2] isa ClassicalOracle
            @test oracle.dcglp isa Model
        end
    end

    @testset "DCGLP Model Structure" begin
        @testset "DCGLP variable existence" begin
            # Build model
            model = Model()
            @variable(model, x[1:3], Bin)
            @variable(model, y[1:3, 1:2] >= 0)
            @objective(model, Min, 10*x[1] + 15*x[2] + 12*x[3] + sum(y))
            @constraint(model, sum(x) >= 1)
            @constraint(model, demand[j in 1:2], sum(y[:,j]) == 1)
            @constraint(model, capacity[i in 1:3, j in 1:2], y[i,j] <= x[i])

            # Decompose
            data, master, oracle = auto_decompose(model, :disjunctive)
            dcglp = oracle.dcglp

            # Verify all DCGLP variables exist
            @test haskey(object_dictionary(dcglp), :tau)
            @test haskey(object_dictionary(dcglp), :omega_0)
            @test haskey(object_dictionary(dcglp), :omega_x)
            @test haskey(object_dictionary(dcglp), :omega_t)
            @test haskey(object_dictionary(dcglp), :sx)
            @test haskey(object_dictionary(dcglp), :st)
        end

        @testset "DCGLP variable dimensions" begin
            # Build model with known dimensions
            model = Model()
            @variable(model, x[1:3], Bin)
            @variable(model, y[1:3, 1:2] >= 0)
            @objective(model, Min, 10*x[1] + 15*x[2] + 12*x[3] + sum(y))
            @constraint(model, sum(x) >= 1)
            @constraint(model, demand[j in 1:2], sum(y[:,j]) == 1)
            @constraint(model, capacity[i in 1:3, j in 1:2], y[i,j] <= x[i])

            # Decompose (dim_x=3, dim_t=1)
            data, master, oracle = auto_decompose(model, :disjunctive)
            dcglp = oracle.dcglp

            # Verify omega_x dimensions [1:2, 1:dim_x]
            omega_x = dcglp[:omega_x]
            @test size(omega_x) == (2, 3)

            # Verify omega_t dimensions [1:2, 1:dim_t]
            omega_t = dcglp[:omega_t]
            @test size(omega_t) == (2, 1)

            # Verify omega_0 dimensions [1:2]
            omega_0 = dcglp[:omega_0]
            @test length(omega_0) == 2

            # Verify sx dimensions [1:dim_x]
            sx = dcglp[:sx]
            @test length(sx) == 3

            # Verify st dimensions [1:dim_t]
            st = dcglp[:st]
            @test length(st) == 1
        end
    end

    @testset "DCGLP Master Constraints" begin
        @testset "Master constraints added to DCGLP" begin
            # Build model with explicit master constraint
            model = Model()
            @variable(model, x[1:3], Bin)
            @variable(model, y[1:3, 1:2] >= 0)
            @objective(model, Min, 10*x[1] + 15*x[2] + 12*x[3] + sum(y))
            @constraint(model, sum(x) >= 1)  # This is a master constraint
            @constraint(model, demand[j in 1:2], sum(y[:,j]) == 1)
            @constraint(model, capacity[i in 1:3, j in 1:2], y[i,j] <= x[i])

            # Classify variables and constraints
            master_vars, oracle_vars = BendersBase.classify_variables_for_benders(model, Symbol[])
            master_constraints, oracle_constraints, coupling_constraints =
                BendersBase.partition_constraints_for_benders(model, master_vars, oracle_vars)

            # Verify we have master constraints
            @test length(master_constraints) > 0

            # Decompose
            data, master, oracle = auto_decompose(model, :disjunctive)
            dcglp = oracle.dcglp

            # Verify DCGLP has constraints (base constraints + master constraints in disjunctive form)
            # Base DCGLP constraints include: omega_t bounds, coneta, condelta, conineq, con0, conx, cont
            # Each master constraint adds 2 constraints (one for kappa, one for nu)
            num_dcglp_constraints = num_constraints(dcglp, count_variable_in_set_constraints=false)
            @test num_dcglp_constraints > 0  # Should have constraints
        end
    end

    @testset "Edge Cases" begin
        @testset "No master constraints" begin
            # Build model where all constraints involve both x and y (coupling constraints only)
            model = Model()
            @variable(model, x[1:2], Bin)
            @variable(model, y[1:2] >= 0)
            @objective(model, Min, sum(x) + sum(y))
            @constraint(model, coupling[i in 1:2], y[i] <= x[i])

            # Decompose should still work
            data, master, oracle = auto_decompose(model, :disjunctive)

            # Verify structure is correct
            @test oracle isa DisjunctiveOracle
            @test length(oracle.typical_oracles) == 2
            @test oracle.dcglp isa Model
        end

        @testset "Custom parameters" begin
            # Build model
            model = Model()
            @variable(model, x[1:2], Bin)
            @variable(model, y[1:2] >= 0)
            @objective(model, Min, sum(x) + sum(y))
            @constraint(model, sum(x) >= 1)
            @constraint(model, linking[i in 1:2], y[i] <= x[i])

            # Create custom parameters
            custom_oracle_param = DisjunctiveOracleParam(
                strengthened = false,
                lift = true
            )
            custom_dcglp_param = DcglpParam(
                time_limit = 500.0,
                verbose = false
            )

            # Decompose with custom parameters
            data, master, oracle = auto_decompose(model, :disjunctive;
                oracle_param = custom_oracle_param,
                dcglp_param = custom_dcglp_param
            )

            # Verify parameters are correctly passed
            @test oracle.oracle_param.strengthened == false
            @test oracle.oracle_param.lift == true
            @test oracle.param.time_limit == 500.0
            @test oracle.param.verbose == false
        end
    end

    @testset "Error Handling" begin
        @testset "Empty model" begin
            # Create empty model with no variables
            model = Model()

            # Should throw ArgumentError
            @test_throws ArgumentError auto_decompose(model, :disjunctive)
        end

        @testset "Wrong oracle_type" begin
            # Create a normal model
            model = Model()
            @variable(model, x[1:2], Bin)
            @variable(model, y[1:2] >= 0)
            @objective(model, Min, sum(x) + sum(y))
            @constraint(model, sum(x) >= 1)

            # Call with :classical should throw error (this method only supports :disjunctive)
            @test_throws ArgumentError auto_decompose(model, :classical)
        end

        @testset "Unknown oracle_type" begin
            # Create a normal model
            model = Model()
            @variable(model, x[1:2], Bin)
            @variable(model, y[1:2] >= 0)
            @objective(model, Min, sum(x) + sum(y))
            @constraint(model, sum(x) >= 1)

            # Call with unknown oracle_type
            @test_throws ArgumentError auto_decompose(model, :unknown)
        end
    end

end

# Old demo code has been replaced with formal tests above
# Previous version used println statements for manual verification
# New version uses @test macros for automated testing
