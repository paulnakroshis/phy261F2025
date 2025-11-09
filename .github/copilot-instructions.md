# AI Agent Instructions for PHY261 Computational Physics

## Project Overview
This repository contains computational physics simulations and exercises for PHY261 at the University of Southern Maine. The codebase primarily uses Julia for numerical simulations of physical systems.

## Key Components and Patterns

### Data Structures
- Physical bodies are represented using immutable `struct Body` with mass, position, and velocity fields
- System states use mutable structs (e.g., `TwoBodySystem`) to track multiple bodies and simulation time
- Vectors are used for position and velocity in 3D space (implemented via `Vector{Float64}`)

### Numerical Methods
- Velocity Verlet integration is used for time evolution of physical systems
- Energy conservation is used as a key validation metric for simulations
- Standard physics constants (like gravitational constant G) are declared as global constants

### Code Organization
- Each simulation is typically contained in a single notebook
- Simulations follow a pattern of:
  1. Data structure definitions
  2. Force/acceleration calculations
  3. Integration methods
  4. Energy conservation checks
  5. Simulation runner with configurable parameters

### Julia-Specific Patterns
- Use `!` suffix for functions that modify their arguments (e.g., `step!`)
- Leverage multiple dispatch for physics calculations
- Use type annotations for performance-critical code
- Employ the LinearAlgebra package for vector operations

## Development Workflow

### Running Simulations
1. Open the notebook in VS Code with Julia extension
2. Execute cells sequentially to load definitions
3. Run simulation functions to see results
4. Check energy conservation metrics for accuracy

### Common Parameters to Adjust
- Time step (`dt`): Controls simulation accuracy vs. speed
- Number of steps (`n_steps`): Controls simulation duration
- Initial conditions: Mass, position, and velocity of bodies
- Output frequency: How often to print/save results

## Best Practices
- Verify energy conservation to validate simulation accuracy
- Use SI units consistently throughout simulations
- Include physical constants with full precision
- Document units in comments for all physical quantities
- Print progress and validation metrics during long simulations

## Limitations and Considerations
- Simulations assume point masses
- Current implementations prioritize clarity over performance
- No relativistic effects are considered
- Simple text-based output (no visualization yet)