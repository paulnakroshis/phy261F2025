using LinearAlgebra
using Plots

# This code will simulate the Earth-Moon system 
# using the Euler-Cromer method for numerical integration.

# Define a structure to hold the state of a body
mutable struct Body
    name::String
    m::Float64
    r::Vector{Vector{Float64}}
    v::Vector{Vector{Float64}}
end

# Gravitational constant
const G = 6.67430e-11  # m^3 kg^-1 s^-2
const day = 86400.0     # seconds in a day
const year = 365.25 * day  # seconds in a year
const AU = 1.496e11    # meters
const Earth_mass = 5.972e24  # kg
const Moon_mass = 7.348e22   # kg
const Earth_radius = 6.371e6  # meters
const Moon_radius = 1.737e6    # meters
const Earth_Moon_distance = 3.844e8  # meters
const Moon_orbital_velocity = 1022.0  # m/s
const total_time = 90.0 * day  # total simulation time
const Δt = 60.0  # time step in seconds
const num_steps = Int(total_time / Δt)
const output_interval = Int(day / Δt)  # output every day
const num_bodies = 2
const body_names = ["Earth", "Moon"]
const body_masses = [Earth_mass, Moon_mass]
const body_radii = [Earth_radius, Moon_radius]
const body_initial_positions = [ [0.0, 0.0, 0.0], [Earth_Moon_distance, 0.0, 0.0] ]
const body_initial_velocities = [ [0.0, 0.0, 0.0], [0.0, Moon_orbital_velocity, 0.0] ]

# Initialize bodies
bodies = Body[]
for i in 1:num_bodies
    push!(bodies, 
          Body(body_names[i], 
          body_masses[i], 
          [body_initial_positions[i]], 
          [body_initial_velocities[i]])
        )
end

# Function to compute the change in velocities for each body due to gravitational forces
function Δv(b1::Body, b2::Body)    
    r_vec = b2.r[end] .- b1.r[end]
    r_mag = norm(r_vec)
    force = G * b1.m * b2.m / r_mag^2
    Δv1 = force * Δt *(r_vec / r_mag) / b1.m
    Δv2 = -force * Δt *(r_vec / r_mag) / b2.m
    return Δv1, Δv2
end

# now use the Euler-Cromer method to update positions and velocities
function EC_Method!(bodies::Vector{Body})    
    for step in 2:num_steps
        # Compute all pairwise gravitational interactions
        for i in 1:length(bodies)-1
            for j in i+1:length(bodies)
                Δv1, Δv2 = Δv(bodies[i], bodies[j])
                push!(bodies[i].v, bodies[i].v[end] .+ Δv1)
                push!(bodies[j].v, bodies[j].v[end] .+ Δv2)
            end
        end
        # Update positions based on new velocities
        for body in bodies
            push!(body.r, body.r[end] .+ body.v[end] .* Δt)
        end
        
        # Output positions at specified intervals
        if step % output_interval == 0
            println("Day $(step * Δt / day):")
            for body in bodies
                println("$(body.name): Position = $(body.r[end]), Velocity = $(body.v[end])")
            end
            println()
        end
    end
end

# Run the simulation and plot the results
EC_Method!(bodies)
# Plotting the trajectories
p = plot()
for body in bodies
    xs = [pos[1] for pos in body.r]
    ys = [pos[2] for pos in body.r]
    plot!(xs, ys, label=body.name)
end
xlabel!("x (m)")
ylabel!("y (m)")
title!("Earth-Moon System Trajectories")
savefig("earth_moon_trajectory.png")
println("Plot saved to earth_moon_trajectory.png")
