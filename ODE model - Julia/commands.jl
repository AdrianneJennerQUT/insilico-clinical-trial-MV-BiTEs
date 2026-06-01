using DifferentialEquations
using Plots

# Parameters
r = 0.2          # growth rate
C₀ = 10.0        # initial population

# ODE definition
function growth!(dC, C, p, t)
    r = p
    dC[1] = r * C[1]
end

# Initial condition and time span
u0 = [C₀]
tspan = (0.0, 20.0)

# Define and solve the problem
prob = ODEProblem(growth!, u0, tspan, r)
sol = solve(prob)

# Plot numerical solution
plot(sol,
     xlabel = "Time",
     ylabel = "C(t)",
     label = "Numerical solution",
     linewidth = 2)

# Plot analytical solution for comparison
t = range(0, 20, length=200)
plot!(t, C₀ .* exp.(r .* t),
      linestyle = :dash,
      label = "Exact solution")

display(plot!)