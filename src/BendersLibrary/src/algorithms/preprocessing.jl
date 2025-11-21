export DisjunctiveRootNodePreprocessing

"""
    DisjunctiveRootNodePreprocessing <: AbstractRootNodePreprocessing

Root node preprocessing using a disjunctive oracle for generating stronger initial cuts.

# Fields
- `classical_oracle::AbstractTypicalOracle`: Oracle for the first phase (classical preprocessing)
- `disjunctive_oracle::AbstractDisjunctiveOracle`: Oracle for the second phase (disjunctive preprocessing)
- `seq_type::Type{<:AbstractBendersSeq}`: Type of sequential Benders algorithm to use
- `params::AbstractBendersSeqParam`: Parameters for the sequential algorithm

# Constructor
```julia
DisjunctiveRootNodePreprocessing(
    classical_oracle::AbstractTypicalOracle,
    disjunctive_oracle::AbstractDisjunctiveOracle;
    seq_type::Type{<:AbstractBendersSeq} = BendersSeq,
    params::AbstractBendersSeqParam = BendersSeqParam()
)
```

# Examples
```julia
# Create classical oracle for first phase
classical_oracle = ClassicalOracle(data)

# Create disjunctive oracle for second phase
disj_oracle = DisjunctiveOracle(data, [oracle_kappa, oracle_nu])

# Create two-phase preprocessing
preprocessing = DisjunctiveRootNodePreprocessing(classical_oracle, disj_oracle)

# Use with BendersBnB
algorithm = BendersBnB(data, master, preprocessing, lazy_callback, user_callback)
```

See also: [`RootNodePreprocessing`](@ref), [`DisjunctiveOracle`](@ref)
"""
mutable struct DisjunctiveRootNodePreprocessing <: AbstractRootNodePreprocessing
    classical_oracle::AbstractTypicalOracle
    disjunctive_oracle::AbstractDisjunctiveOracle
    seq_type::Type{<:AbstractBendersSeq}
    params::AbstractBendersSeqParam

    function DisjunctiveRootNodePreprocessing(
        classical_oracle::AbstractTypicalOracle,
        disjunctive_oracle::AbstractDisjunctiveOracle;
        seq_type::Type{<:AbstractBendersSeq} = BendersSeq,
        params::AbstractBendersSeqParam = BendersSeqParam()
    )
        new(classical_oracle, disjunctive_oracle, seq_type, params)
    end
end

"""
    root_node_processing!(data::Data, master::AbstractMaster, preprocessing::DisjunctiveRootNodePreprocessing) -> Float64

Perform two-phase root node preprocessing.

This method relaxes integrality constraints once and executes two sequential phases:
1. Classical oracle preprocessing to generate initial cuts
2. Disjunctive oracle preprocessing to generate strengthened cuts

The time limit for the second phase is automatically adjusted based on the time
consumed by the first phase.

# Arguments
- `data::Data`: Problem data
- `master::AbstractMaster`: Master problem
- `preprocessing::DisjunctiveRootNodePreprocessing`: Two-phase preprocessing configuration

# Returns
- `Float64`: Total time spent on both preprocessing phases
"""
function root_node_processing!(data::Data, master::AbstractMaster, preprocessing::DisjunctiveRootNodePreprocessing)
    total_time = 0.0

    # Relax integrality once for both phases
    undo = relax_integrality(master.model)

    # Phase 1: Classical oracle preprocessing
    classical_time = @elapsed begin
        classical_param = deepcopy(preprocessing.params)
        classical_seq = preprocessing.seq_type(data, master, preprocessing.classical_oracle; param=classical_param)
        solve!(classical_seq)
    end
    total_time += classical_time

    # Phase 2: Disjunctive oracle preprocessing with remaining time
    remaining_param = deepcopy(preprocessing.params)
    remaining_param.time_limit -= classical_time

    disjunctive_time = @elapsed begin
        disjunctive_seq = preprocessing.seq_type(data, master, preprocessing.disjunctive_oracle; param=remaining_param)
        solve!(disjunctive_seq)
    end
    total_time += disjunctive_time

    # Restore integrality
    undo()

    return total_time
end
