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
            @test length(oracle.fixed_x_constraints) == 3  # One fixing constraint per master variable
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

    @testset "Explicit decision variables with fractional master variables" begin
        @testset "Basic decision_vars functionality" begin
            # Build model with only continuous variables
            model = Model()
            @variable(model, capacity[1:2] >= 0)
            @variable(model, y[1:2, 1:3] >= 0)
            @objective(model, Min, sum(capacity) + sum(y))
            @constraint(model, sum(capacity) >= 10)
            @constraint(model, demand[j in 1:3], sum(y[:,j]) == 1)

            # Classify with capacity explicitly as decision variable
            master_vars, oracle_vars = BendersBase.classify_variables_for_benders(model, [Symbol("capacity[1]"), Symbol("capacity[2]")])

            # Verify capacity variables are in master
            @test length(master_vars) == 2
            @test length(oracle_vars) == 6

            # Verify the correct variables are classified
            master_var_names = Set(name(v) for v in master_vars)
            @test "capacity[1]" in master_var_names
            @test "capacity[2]" in master_var_names
        end

        @testset "Mixed automatic and explicit classification" begin
            # Build model with binary and continuous variables
            model = Model()
            @variable(model, x[1:2], Bin)
            @variable(model, capacity[1:2] >= 0)
            @variable(model, y[1:2, 1:3] >= 0)
            @objective(model, Min, sum(x) + sum(capacity) + sum(y))
            @constraint(model, sum(x) >= 1)
            @constraint(model, sum(capacity) >= 10)
            @constraint(model, demand[j in 1:3], sum(y[:,j]) == 1)

            # Classify with capacity explicitly as decision variable
            master_vars, oracle_vars = BendersBase.classify_variables_for_benders(model, [Symbol("capacity[1]"), Symbol("capacity[2]")])

            # Verify both x (automatic) and capacity (explicit) are in master
            @test length(master_vars) == 4  # 2 x + 2 capacity
            @test length(oracle_vars) == 6  # 6 y

            # Verify the correct variables are classified
            master_var_names = Set(name(v) for v in master_vars)
            @test "x[1]" in master_var_names
            @test "x[2]" in master_var_names
            @test "capacity[1]" in master_var_names
            @test "capacity[2]" in master_var_names
        end

        @testset "Decision_vars with constraint partitioning" begin
            # Build model
            model = Model()
            @variable(model, x[1:2], Bin)
            @variable(model, capacity[1:2] >= 0)
            @variable(model, y[1:2, 1:3] >= 0)
            @objective(model, Min, sum(x) + sum(capacity) + sum(y))
            @constraint(model, sum(x) >= 1)
            @constraint(model, sum(capacity) >= 10)
            @constraint(model, demand[j in 1:3], sum(y[:,j]) == 1)

            # Classify and partition
            master_vars, oracle_vars = BendersBase.classify_variables_for_benders(model, [Symbol("capacity[1]"), Symbol("capacity[2]")])
            master_constraints, oracle_constraints, coupling_constraints =
                BendersBase.partition_constraints_for_benders(model, master_vars, oracle_vars)

            # Verify constraint partitioning
            @test length(master_constraints) == 2  # sum(x) >= 1 and sum(capacity) >= 10
            @test length(oracle_constraints) == 3  # 3 demand constraints
            @test length(coupling_constraints) == 0  # No coupling constraints

            # Verify total count
            total_constraints = length(master_constraints) + length(oracle_constraints) + length(coupling_constraints)
            @test total_constraints == num_constraints(model, count_variable_in_set_constraints=false)
        end

        @testset "Decision_vars with objective decomposition" begin
            # Build model
            model = Model()
            @variable(model, x[1:2], Bin)
            @variable(model, capacity[1:2] >= 0)
            @variable(model, y[1:2, 1:3] >= 0)
            @objective(model, Min, sum(x) + sum(capacity) + sum(y))

            # Classify and decompose
            master_vars, oracle_vars = BendersBase.classify_variables_for_benders(model, [Symbol("capacity[1]"), Symbol("capacity[2]")])
            master_objective, oracle_objective =
                BendersBase.decompose_objective_for_benders(model, master_vars, oracle_vars)

            # Verify master objective contains x and capacity
            @test length(master_objective.terms) == 4  # 2 x + 2 capacity
            @test length(oracle_objective.terms) == 6  # 6 y

            # Verify capacity variables are in master objective
            capacity_in_master = false
            for (var, coeff) in master_objective.terms
                if startswith(name(var), "capacity")
                    capacity_in_master = true
                    @test coeff == 1.0
                end
            end
            @test capacity_in_master
        end

        @testset "Full auto_decompose with decision_vars" begin
            # Build model with coupling constraints
            model = Model()
            @variable(model, x[1:2], Bin)
            @variable(model, 0 <= capacity[1:2] <= 100)
            @variable(model, y[1:2, 1:3] >= 0)
            @objective(model, Min, 10*sum(x) + 5*sum(capacity) + sum(y))
            @constraint(model, sum(x) >= 1)
            @constraint(model, sum(capacity) >= 10)
            @constraint(model, demand[j in 1:3], sum(y[:,j]) == 1)
            @constraint(model, link[i in 1:2, j in 1:3], y[i,j] <= capacity[i])

            # Decompose with decision_vars
            data, master, oracle = auto_decompose(model; decision_vars=[Symbol("capacity[1]"), Symbol("capacity[2]")])

            # Verify Data structure
            @test data.dim_x == 4  # 2 x + 2 capacity
            @test data.dim_t == 1
            @test length(data.c_x) == 4

            # Verify coefficient values (x variables have coeff 10, capacity have coeff 5)
            # Note: The order depends on the order returned by all_variables(model)
            # We verify that we have two 10.0s and two 5.0s
            coeff_counts = Dict(10.0 => 0, 5.0 => 0)
            for c in data.c_x
                if c == 10.0
                    coeff_counts[10.0] += 1
                elseif c == 5.0
                    coeff_counts[5.0] += 1
                end
            end
            @test coeff_counts[10.0] == 2
            @test coeff_counts[5.0] == 2

            # Verify Master model structure
            @test num_variables(master.model) == 5  # 4 x/capacity + 1 t
            @test length(master.x) == 4
            @test length(master.t) == 1

            # Verify Oracle model structure
            @test length(oracle.fixed_x_constraints) == 4  # 4 master variables
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

        @testset "Quadratic constraint rejection" begin
            # Create model with quadratic constraint
            model = Model()
            @variable(model, x[1:2], Bin)
            @variable(model, y[1:2] >= 0)
            @objective(model, Min, sum(x) + sum(y))
            @constraint(model, x[1] + y[1] >= 1)
            @constraint(model, y[1]^2 + y[2]^2 <= 1)  # Quadratic constraint

            # Should throw ArgumentError
            @test_throws ArgumentError auto_decompose(model)
        end

        @testset "Nonlinear constraint rejection" begin
            # Create model with nonlinear constraint
            model = Model()
            @variable(model, x[1:2], Bin)
            @variable(model, y[1:2] >= 0)
            @objective(model, Min, sum(x) + sum(y))
            @constraint(model, x[1] + y[1] >= 1)
            @constraint(model, y[1] * y[2] <= 1)  # Nonlinear constraint

            # Should throw ArgumentError
            @test_throws ArgumentError auto_decompose(model)
        end
    end

end


