# Architecture Design Proposal: Modular Subpackage Structure for BendersDecomposition.jl

---

## Executive Summary

This proposal presents a comprehensive architectural comparison between two design approaches for migrating problem implementations from the `example/` directory into the core package structure. After thorough analysis, **we recommend the modular subpackage architecture (Approach 2)** featuring `BendersBase.jl` and `BendersLibrary.jl` as independent submodules. This design offers superior modularity, maintainability, dependency isolation, and aligns with Julia ecosystem best practices for scientific computing packages.

### Key Findings

| Criterion | Approach 1 (Monolithic) | Approach 2 (Modular) | Winner |
|-----------|------------------------|----------------------|---------|
| **Modularity** | Low - tight coupling | High - clear separation | ✓ Approach 2 |
| **Dependency Management** | Mixed - all dependencies loaded | Isolated - minimal core deps | ✓ Approach 2 |
| **Testing Isolation** | Difficult - interdependencies | Clean - independent test suites | ✓ Approach 2 |
| **Extensibility** | Moderate - namespace conflicts | Excellent - plugin architecture | ✓ Approach 2 |
| **Maintenance Burden** | High - monolithic codebase | Low - separate concerns | ✓ Approach 2 |
| **Initial Complexity** | Low - simpler structure | Moderate - requires coordination | Approach 1 |

**Recommendation**: Adopt Approach 2 (Modular Subpackage Architecture)

---

## 1. Background and Context

### 1.1 Problem Statement

The original `asu-opt-lab/BendersDecomposition` repository contains problem implementations in an `example/` directory, including:

- **UFLP** (Uncapacitated Facility Location Problem)
- **CFLP** (Capacitated Facility Location Problem)
- **SCFLP** (Stochastic Capacitated Facility Location Problem)
- **MCNDP** (Multi-Commodity Network Design Problem)
- **SNIP** (Stochastic Network Interdiction Problem)

These implementations currently exist as standalone examples rather than first-class package components. We need to migrate them into the formal package structure while maintaining code quality, testability, and extensibility.

### 1.2 Original Repository Structure

```
asu-opt-lab/BendersDecomposition/
├── src/
│   ├── BendersDecomposition.jl     # Main module
│   ├── types.jl                     # Core abstract types
│   ├── algorithms/                  # Core algorithms
│   │   ├── BendersSeq.jl
│   │   ├── BendersBnB.jl
│   │   ├── Dcglp.jl
│   │   └── callback/
│   ├── modules/                     # Core components
│   │   ├── master.jl
│   │   ├── oracle.jl
│   │   └── oracleTypical*.jl
│   └── utils/                       # Utility functions
│       ├── utilsSolver.jl
│       ├── utilsCuts.jl
│       └── utilsLoop*.jl
└── example/                         # Problem implementations
    ├── UFLP/
    ├── CFLP/
    ├── SCFLP/
    ├── MCNDP/
    └── SNIP/
```

---

## 2. Proposed Architectural Approaches

### 2.1 Approach 1: Monolithic Package with `problems/` Directory

**Description**: Directly migrate `example/` into `src/problems/` as a subdirectory within a single package namespace.

```julia
src/
├── BendersDecomposition.jl         # Main module with namespace exports
├── types.jl                         # Core abstract types
│
├── algorithms/                      # Core algorithms (unchanged)
│   ├── algorithms.jl
│   ├── BendersSeq.jl
│   ├── BendersSeqInOut.jl
│   ├── BendersBnB.jl
│   ├── Dcglp.jl
│   └── callback/
│       ├── BendersCallback.jl
│       └── LazyCallback.jl
│
├── modules/                         # Core components (unchanged)
│   ├── modules.jl
│   ├── master.jl
│   ├── oracle.jl
│   └── oracleTypical*.jl
│
├── utils/                           # Utilities (unchanged)
│   ├── utils.jl
│   ├── utilsSolver.jl
│   ├── utilsCuts.jl
│   └── utilsLoop*.jl
│
└── problems/                        # NEW: Problem implementations
    ├── UFLP/
    │   ├── UFLP.jl                 # Submodule entry point
    │   ├── data.jl                 # UFLPData + readers
    │   ├── oracles.jl              # UFLKnapsackOracle
    │   └── models.jl               # update_model! overloads
    ├── CFLP/
    │   ├── CFLP.jl
    │   ├── data.jl
    │   ├── oracles.jl
    │   └── models.jl
    ├── SCFLP/
    ├── SNIP/
    └── MCNDP/
```

**Characteristics**:
- Single `Project.toml` with unified dependencies
- All code under `BendersDecomposition` namespace
- Problem implementations as nested modules
- Shared utility functions across all components

---

### 2.2 Approach 2: Modular Subpackage Architecture (Recommended)

**Description**: Separate the package into two independent submodules with distinct responsibilities and dependency chains.

```julia
BendersDecomposition.jl/             # Meta-package
├── Project.toml                     # Umbrella package dependencies
├── src/
│   └── BendersDecomposition.jl     # Re-exports BendersBase + BendersLibrary
│
├── src/BendersBase/                 # SUBPACKAGE 1: Core Framework
│   ├── Project.toml                 # Minimal deps: JuMP, MOI, MacroTools
│   ├── src/
│   │   ├── BendersBase.jl          # Core module
│   │   ├── types.jl                 # Abstract type hierarchy
│   │   ├── modeling/
│   │   │   ├── modeling.jl
│   │   │   ├── macros.jl            # @benders_decomposition DSL
│   │   │   ├── stochastic_macros.jl
│   │   │   └── auto_decomposition.jl
│   │   ├── solution_procedure/
│   │   │   ├── sequential.jl
│   │   │   ├── callbacks.jl
│   │   │   ├── state.jl
│   │   │   └── io.jl
│   │   ├── cut_strategy/
│   │   │   ├── cuts.jl              # ClassicalCut, ParetoCut, UnifiedCut
│   │   │   └── stochastic_cut.jl
│   │   └── config.jl
│   └── test/                        # Comprehensive unit tests
│
└── src/BendersLibrary/              # SUBPACKAGE 2: Problem Implementations
    ├── Project.toml                 # Deps: BendersBase, Distributions, JSON
    ├── src/
    │   ├── BendersLibrary.jl        # Problem library module
    │   ├── types.jl                 # Concrete data types (CFLPData, etc.)
    │   ├── modeling/
    │   │   ├── modeling.jl
    │   │   ├── UFLP.jl              
    │   │   ├── CFLP.jl
    │   │   ├── SCFLP.jl
    │   │   ├── MCNDP.jl             
    │   │   └── SNIP.jl              
    │   ├── cut/
    │   │   ├── knapsack_cut.jl      # KnapsackCut strategy
    │   │   └── fs_knapsack_cut.jl   # FatKnapsackCut, SlimKnapsackCut
    │   ├── solution_procedure/
    │   │   └── StabilizedSequential.jl
    │   └── utilities/
    │       ├── data_reader.jl
    │       └── data_generator.jl
    └── test/                        # Problem-specific integration tests
```

**Characteristics**:
- Two independent `Project.toml` files with separate dependency chains
- `BendersBase`: Minimal dependencies (JuMP + MOI only)
- `BendersLibrary`: Additional dependencies (Distributions, JSON, etc.)
- Clear separation of concerns: framework vs. applications
- Independent version control and release cycles

---

## 3. Detailed Comparative Analysis

### 3.1 Modularity and Separation of Concerns

#### Approach 1: Monolithic Structure
- **Coupling**: High coupling between core algorithms and problem implementations
- **Namespace**: Single flat namespace `BendersDecomposition.UFLP`, `BendersDecomposition.CFLP`
- **Abstraction Boundaries**: Blurred boundaries between framework and applications
- **Code Organization**: All code lives in one package, requiring careful manual namespace management

**Assessment**: ⚠️ Risk of tight coupling and namespace pollution as the codebase grows.

#### Approach 2: Modular Subpackages
- **Coupling**: Loose coupling via well-defined abstract interfaces
- **Namespace**: Clear hierarchical namespaces: `BendersBase.*` vs. `BendersLibrary.*`
- **Abstraction Boundaries**: Strict separation enforced by package boundaries
- **Code Organization**: Natural division of responsibilities

**Assessment**: ✅ Superior modularity with clear architectural boundaries.

---

### 3.2 Dependency Management

#### Approach 1: Unified Dependencies
```toml
[deps]
JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
MacroTools = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"      # Only needed for SCFLP
JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"               # Only needed for data I/O
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"             # Only needed for MCNDP
# ... all dependencies mixed together
```

**Problems**:
- Users who only need the core framework must install problem-specific dependencies
- Dependency bloat for minimal use cases
- Harder to track which dependencies are actually required for core functionality
- Increased installation time and potential version conflicts

#### Approach 2: Isolated Dependencies

**BendersBase/Project.toml** (Minimal Core):
```toml
[deps]
JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
MacroTools = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
# Only essential dependencies - no problem-specific packages
```

**BendersLibrary/Project.toml** (Problem-Specific):
```toml
[deps]
BendersBase = "a1b2c3d4-e5f6-7890-1234-567890abcdef"       # Core framework
JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"      # For stochastic problems
JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"               # For data I/O
# Problem-specific dependencies isolated here
```

**Benefits**:
- Users can install only `BendersBase` for custom implementations
- Problem-specific dependencies isolated in `BendersLibrary`
- Cleaner dependency graph and faster precompilation
- Easier maintenance and security updates

**Assessment**: ✅ Approach 2 provides superior dependency isolation and user flexibility.

---

### 3.3 Testing and Quality Assurance

#### Approach 1: Unified Test Suite
```julia
test/
├── runtests.jl                      # All tests in one suite
├── test_algorithms.jl
├── test_master.jl
├── test_oracle.jl
├── test_uflp.jl                     # Mixed with core tests
├── test_cflp.jl
└── test_scflp.jl
```

**Challenges**:
- Core framework tests mixed with problem-specific tests
- Difficult to run only framework tests without problem dependencies
- Test failures in problem implementations can block core development
- Harder to maintain independent CI/CD pipelines

#### Approach 2: Isolated Test Suites

**BendersBase/test/**:
```julia
test/
├── runtests.jl                      # Framework tests only
├── test_types.jl
├── modeling/
│   ├── test_modeling.jl
│   ├── test_macros.jl
│   └── test_auto_decomposition.jl
├── solution_procedure/
│   ├── test_sequential.jl
│   ├── test_callbacks.jl
│   └── test_state.jl
└── cut_strategy/
    ├── test_cuts.jl
    └── test_stochastic_cut.jl
```

**BendersLibrary/test/**:
```julia
test/
├── runtests.jl                      # Problem-specific integration tests
├── test_uflp.jl
├── test_cflp.jl
├── test_scflp.jl
├── test_knapsack_cuts.jl
└── test_data_readers.jl
```

**Benefits**:
- Independent test execution: `julia --project=src/BendersBase -e "using Pkg; Pkg.test()"`
- Framework tests run without problem-specific dependencies
- Parallel CI/CD: test both packages simultaneously
- Clear test ownership and responsibility

**Assessment**: ✅ Approach 2 enables cleaner testing practices and faster CI/CD.

---

### 3.4 Extensibility and Plugin Architecture

#### Approach 1: Extension with Heavy Dependencies
```julia
# User creates external package: MyBendersProblem.jl
using BendersDecomposition  # Must import entire package

struct MyCustomData <: BendersDecomposition.AbstractData
    # Custom problem data
end

struct MyCustomOracle <: BendersDecomposition.BaseOracle
    # Custom oracle implementation
end

# User's Project.toml:
[deps]
BendersDecomposition = "..."  # Brings ALL dependencies
# ↓ Forced to install even if not needed:
# - Distributions.jl (only for SCFLP)
# - JSON.jl (only for data readers)
# - Graphs.jl (if MCNDP added later)
```

**Limitations**:
- **Dependency Bloat**: External packages must install ALL problem-specific dependencies, even unused ones
- **Tight Version Coupling**: User packages tied to monolithic version updates (if CFLP changes, all extensions affected)
- **Namespace Pollution**: All problem implementations visible in user's namespace
- **Breaking Changes Risk**: Library additions/changes can break external extensions

**Technical Note**: Extension is possible via Julia's multiple dispatch, but the dependency burden makes it impractical for lightweight custom implementations.

#### Approach 2: Lightweight Plugin Architecture
```julia
# User creates independent package: MyBendersProblem.jl
using BendersBase  # Import only the framework

struct MyCustomData <: AbstractData
    # Custom problem data
end

struct MyCustomOracle <: BaseOracle
    # Custom oracle implementation
end

# Implement required interfaces:
BendersBase.generate_cuts(oracle::MyCustomOracle, ...) = ...
BendersBase.solve!(problem::MyCustomProblem, ...) = ...
```

**Advantages**:
- **Minimal Dependencies**: Only `BendersBase` (JuMP + MOI) required - typically 10-20% of monolithic package size
- **Version Independence**: Framework stability independent of problem library evolution
- **Clean Namespace**: Only import what you need (`using BendersBase` vs. `using BendersDecomposition`)
- **Publish Independently**: Research packages can be registered separately on Julia General Registry
- **True Plugin Architecture**: Framework evolution doesn't break external extensions


**Assessment**: ✅ Approach 2 dramatically reduces barrier to entry for external extensions.

---

### 3.5 Solver Compatibility and Defaults

**Critical Issue**: Default solver selection

#### Approach 1: Risk of Commercial Solver Lock-in
In a monolithic structure, problem implementations may inadvertently use commercial solvers (CPLEX, Gurobi) as defaults:

```julia
# Bad practice in examples
function solve_uflp(data)
    model = Model(CPLEX.Optimizer)  # ❌ Commercial solver as default
    # ...
end
```

**Problems**:
- New users cannot run examples without commercial licenses
- Academic users face barriers to entry
- Reduces package accessibility and adoption

#### Approach 2: Clean Solver Abstraction

**BendersBase.jl**: Solver-agnostic framework
```julia
# Framework never specifies solvers
function solve!(problem::AbstractBendersProblem, ...)
    # Assumes optimizer already set by user
    optimize!(problem.master.model)
end
```

**BendersLibrary.jl** and **Documentation**: Use open-source defaults
```julia
# Examples and tutorials use HiGHS (open-source)
using HiGHS
set_optimizer(benders_problem.master.model, HiGHS.Optimizer)
set_optimizer(benders_problem.oracle.model, HiGHS.Optimizer)
```

**Test Configuration**:
```toml
# BendersBase/Project.toml
[extras]
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
HiGHS = "87dc4568-4c63-4d18-b0c0-bb2238e4078b"  # Open-source solver

[targets]
test = ["Test", "HiGHS"]
```

**Benefits**:
- Framework works with any JuMP-compatible solver
- Examples run out-of-the-box without licenses
- Users can easily swap solvers: HiGHS → CPLEX → Gurobi
- Broader accessibility for academic and open-source users

**Assessment**: ✅ Approach 2 naturally enforces solver-agnostic design.

---

### 3.6 Documentation and User Experience

#### Approach 1: Mixed Documentation
```
docs/
├── manual/
│   ├── core_concepts.md              # Framework docs
│   ├── algorithms.md
│   ├── uflp_guide.md                 # Problem-specific docs mixed in
│   └── cflp_guide.md
└── api/
    └── reference.md                  # All APIs in one namespace
```

**Challenges**:
- Users must navigate mixed framework/problem documentation
- API reference includes everything (overwhelming for new users)
- Unclear learning path: "What do I need vs. what's optional?"

#### Approach 2: Layered Documentation

**BendersBase Documentation**:
```markdown
# Getting Started with BendersBase.jl
Learn the core framework concepts:
- Abstract types and interfaces
- Modeling with macros
- Solution procedures
- Custom cut strategies

# API Reference
Only framework APIs - concise and focused
```

**BendersLibrary Documentation**:
```markdown
# Problem Library Guide
Browse problem implementations:
- UFLP Tutorial (with KnapsackCut)
- CFLP Tutorial (with strengthened cuts)
- SCFLP Tutorial (stochastic programming)

# Adding Your Own Problem
Step-by-step guide to extend the library
```

**User Learning Path**:
1. Start with `BendersBase` docs → understand framework
2. Browse `BendersLibrary` examples → learn from implementations
3. Build custom problem using learned patterns

**Assessment**: ✅ Approach 2 provides clearer learning progression.

---


## 4. Industry Best Practices and Ecosystem Alignment

### 4.1 Julia Package Ecosystem Patterns

**Successful Examples of Modular Subpackage Architecture**:

| Framework Package | Algorithm/Solver Packages | Pattern |
|-------------------|---------------------------|---------|
| `DiffEqBase.jl` | `OrdinaryDiffEq.jl`, `StochasticDiffEq.jl` | Core + Specialized Solvers |
| `AbstractPlotting.jl` | `Makie.jl`, `GLMakie.jl`, `CairoMakie.jl` | Core + Backends |
| `Tables.jl` | `DataFrames.jl`, `Arrow.jl`, `CSV.jl` | Interface + Implementations |
| `JuMP.jl` | `GLPK.jl`, `Gurobi.jl`, `CPLEX.jl` | Modeling + Solvers |

**Common Pattern**:
```julia
# Framework defines interfaces
AbstractFramework.jl  →  Minimal deps, stable API

# Implementations provide concrete types
ImplementationPackage.jl  →  Depends on framework, evolves rapidly
```

**BendersDecomposition.jl follows this proven pattern**:
```julia
BendersBase.jl (Framework)  →  Abstract types, core algorithms
BendersLibrary.jl (Implementations)  →  Concrete problems, specialized cuts
```

### 4.2 SciML Ecosystem Alignment

Our `BendersBase` + `BendersLibrary` structure mirrors the **SciML (Scientific Machine Learning) ecosystem design**:

```julia
# SciML Pattern:
DiffEqBase.jl (interfaces) → OrdinaryDiffEq.jl (solvers)

# Our Pattern:
BendersBase.jl (interfaces) → BendersLibrary.jl (problems)
```

**Benefits of Alignment**:
- Users familiar with SciML will immediately understand our architecture
- Interoperability with SciML packages (e.g., sensitivity analysis, optimization)
- Access to SciML's documentation patterns and tooling

---


## 5. Conclusion

After comprehensive analysis across modularity, dependency management, testing, extensibility, maintenance, and ecosystem alignment, **Approach 2 (Modular Subpackage Architecture) is the clear superior choice**.

### Key Advantages Summary

✅ **Technical Excellence**:
- Clean separation of concerns via package boundaries
- Minimal dependency footprint for core framework
- Independent testing and CI/CD pipelines

✅ **Ecosystem Alignment**:
- Follows proven patterns from DiffEqBase, JuMP, Tables
- Enables third-party extensions via plugin architecture
- Solver-agnostic design avoiding commercial lock-in

✅ **Practical Benefits**:
- Already implemented - minimal migration risk
- Clear maintenance responsibilities
- Independent versioning and release cycles
- Superior documentation structure

✅ **Future-Proof Design**:
- Supports external problem packages
- Stable framework with evolving problem library
- Research-friendly extension mechanism

