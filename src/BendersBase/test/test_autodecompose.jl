using Test
using BendersBase
using JuMP

@testset "Auto Decompose" begin

    @testset "Basic Functionality" begin
        @testset "Simple facility location model" begin
            # Build a simple facility location model
            model = Model()
            @variable(model, x[1:3], Bin)
            @variable(model, y[1:3, 1:2] >= 0)
            @objective(model, Min, 10*x[1] + 15*x[2] + 12*x[3] + sum(y))
            @constraint(model, sum(x) >= 1)
            @constraint(model, demand[j in 1:2], sum(y[:,j]) == 1)
            @constraint(model, capacity[i in 1:3, j in 1:2], y[i,j] <= x[i])

            # Decompose the model
            data, master, oracle = auto_decompose(model)

            # Verify Data structure
            @test data.dim_x == 3
            @test data.dim_t == 1
            @test length(data.c_x) == 3
            @test length(data.c_t) == 1
            @test data.c_t[1] == 1.0

            # Verify Master model structure
            @test num_variables(master.model) == 4  # 3 x variables + 1 t variable
            @test length(master.x) == 3
            @test length(master.t) == 1

            # Verify Oracle model structure
            @test length(oracle.fixing_constraints) == 3  # One fixing constraint per master variable
        end
    end

    @testset "Classification Tests" begin
        @testset "Constraint partitioning" begin
            # Build model
            model = Model()
            @variable(model, x[1:3], Bin)
            @variable(model, y[1:3, 1:2] >= 0)
            @objective(model, Min, 10*x[1] + 15*x[2] + 12*x[3] + sum(y))
            @constraint(model, sum(x) >= 1)
            @constraint(model, demand[j in 1:2], sum(y[:,j]) == 1)
            @constraint(model, capacity[i in 1:3, j in 1:2], y[i,j] <= x[i])

            # Classify variables
            master_vars, oracle_vars = BendersBase.classify_variables_for_benders(model, Symbol[])

            # Partition constraints
            master_constraints, oracle_constraints, coupling_constraints =
                BendersBase.partition_constraints_for_benders(model, master_vars, oracle_vars)

            # Verify constraint partitioning
            @test length(master_constraints) == 1  # sum(x) >= 1
            @test length(oracle_constraints) == 2  # 2 demand constraints
            @test length(coupling_constraints) == 6  # 6 capacity constraints

            # Verify total count
            total_constraints = length(master_constraints) + length(oracle_constraints) + length(coupling_constraints)
            @test total_constraints == num_constraints(model, count_variable_in_set_constraints=false)
        end

        @testset "Objective decomposition" begin
            # Build model
            model = Model()
            @variable(model, x[1:3], Bin)
            @variable(model, y[1:3, 1:2] >= 0)
            @objective(model, Min, 10*x[1] + 15*x[2] + 12*x[3] + sum(y))

            # Classify variables
            master_vars, oracle_vars = BendersBase.classify_variables_for_benders(model, Symbol[])

            # Decompose objective
            master_objective, oracle_objective =
                BendersBase.decompose_objective_for_benders(model, master_vars, oracle_vars)

            # Verify master objective contains x coefficients
            @test haskey(master_objective.terms, master_vars[1])
            @test master_objective.terms[master_vars[1]] == 10.0
            @test master_objective.terms[master_vars[2]] == 15.0
            @test master_objective.terms[master_vars[3]] == 12.0

            # Verify oracle objective contains y coefficients
            @test length(oracle_objective.terms) == 6  # 3*2 y variables
        end
    end

    @testset "Edge Cases" begin
        @testset "Pure binary variable model" begin
            # Build model with only binary variables
            model = Model()
            @variable(model, x[1:3], Bin)
            @objective(model, Min, sum(x))
            @constraint(model, sum(x) >= 1)

            # Classify variables
            master_vars, oracle_vars = BendersBase.classify_variables_for_benders(model, Symbol[])

            # Verify all variables are classified as master variables
            @test length(master_vars) == 3
            @test length(oracle_vars) == 0

            # Decompose should still work
            data, master, oracle = auto_decompose(model)
            @test data.dim_x == 3
        end

        @testset "Oracle-only objective" begin
            # Build model where objective only contains continuous variables
            model = Model()
            @variable(model, x[1:2], Bin)
            @variable(model, y[1:3] >= 0)
            @objective(model, Min, sum(y))  # Only y variables in objective
            @constraint(model, sum(x) >= 1)
            @constraint(model, linking[i in 1:2], y[i] <= x[i])

            # Classify variables
            master_vars, oracle_vars = BendersBase.classify_variables_for_benders(model, Symbol[])

            # Decompose objective
            master_objective, oracle_objective =
                BendersBase.decompose_objective_for_benders(model, master_vars, oracle_vars)

            # Verify master objective has no variable terms (only constant)
            @test length(master_objective.terms) == 0

            # Verify oracle objective has y variable terms
            @test length(oracle_objective.terms) == 3
        end
    end

    @testset "Error Handling" begin
        @testset "Empty model" begin
            # Create empty model with no variables
            model = Model()

            # Should throw ArgumentError
            @test_throws ArgumentError auto_decompose(model)
        end
    end

end


