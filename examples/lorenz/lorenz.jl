# What the examples on the Lorenz system share: the system at σ = 10, ρ = 28, β = 8/3, sampled every
# 0.05 by five steps of 0.01 of the classical Runge–Kutta scheme, its two distances, the invariant of
# its rotation, and the look of the figures. Each example includes this file,
#
#     include(joinpath(@__DIR__, "lorenz.jl"))
#
# The rotation R(x, y, z) = (-x, -y, z) by π about the z axis maps a solution onto a solution. The
# distance up to it is the smaller of ‖b - a‖ and ‖R b - a‖ over ‖a‖, and it returns the rotations r
# of b, 0 or 1, that reach it. The invariant v = (√(x² + y²), z) does not change under R, has the norm
# of the state, and does not stretch distances, by the reverse triangle inequality on (x, y).
import NearRecurrences
using CairoMakie
using Printf


# ---------------------------------------------------------------------------- #
# the system                                                                   #
# ---------------------------------------------------------------------------- #

function lorenz!(dx, x)
    dx[1] = 10*(x[2] - x[1])
    dx[2] = x[1]*(28 - x[3]) - x[2]
    dx[3] = x[1]*x[2] - 8/3*x[3]
    return dx
end

const WORK = ntuple(_ -> zeros(3), 5)

# ---- one step of 0.01 of the classical scheme, in place ----
function step!(x)
    k₁, k₂, k₃, k₄, y = WORK
    lorenz!(k₁, x)
    y .= x .+ 0.005 .* k₁
    lorenz!(k₂, y)
    y .= x .+ 0.005 .* k₂
    lorenz!(k₃, y)
    y .= x .+ 0.01 .* k₃
    lorenz!(k₄, y)
    x .+= (0.01/6) .* (k₁ .+ 2 .* k₂ .+ 2 .* k₃ .+ k₄)
    return x
end

# ---- one sample ahead, five steps ----
function advance!(x)
    for _ in 1:5
        step!(x)
    end
    return x
end

# ---- the points of the trajectory from u, every step of 0.01 over n steps, for the curves of the
#      figures ----
function fine(u, n)
    x      = copy(u)
    points = [copy(x)]
    for _ in 1:n
        step!(x)
        push!(points, copy(x))
    end
    return points
end

# ---- the state after a transient of 20 time units from (1, 1, 1) ----
function start()
    x = [1.0, 1.0, 1.0]
    for _ in 1:400
        advance!(x)
    end
    return x
end

# ---- ‖b - a‖/‖a‖ ----
distance(a, b) = sqrt((b[1] - a[1])^2 + (b[2] - a[2])^2 + (b[3] - a[3])^2)/
                 sqrt(a[1]^2 + a[2]^2 + a[3]^2)

# ---- the least of ‖b - a‖ and ‖R b - a‖ over ‖a‖, and the rotations r of b that reach it ----
function rotated(a, b)
    d₀ = sqrt((b[1] - a[1])^2 + (b[2] - a[2])^2 + (b[3] - a[3])^2)
    d₁ = sqrt((-b[1] - a[1])^2 + (-b[2] - a[2])^2 + (b[3] - a[3])^2)
    na = sqrt(a[1]^2 + a[2]^2 + a[3]^2)
    return d₁ < d₀ ? (d₁/na, 1) : (d₀/na, 0)
end

# ---- (√(x² + y²), z), unchanged by R, of norm ‖a‖ ----
invariant(a) = [hypot(a[1], a[2]), a[3]]


# ---------------------------------------------------------------------------- #
# the look                                                                     #
# ---------------------------------------------------------------------------- #

const BLUE   = Makie.to_color("#1976D2")
const RED    = Makie.to_color("#E53935")
const PURPLE = Makie.to_color("#5E35B1")
const GREY   = Makie.to_color("#757575")

set_theme!(theme_latexfonts())
update_theme!(fontsize = 17,
              Axis     = (topspinevisible   = false,
                          rightspinevisible = false,
                          xgridvisible      = false,
                          ygridvisible      = false,
                          titlealign        = :left,
                          titlefont         = :regular,
                          titlesize         = 17),
              Legend   = (framevisible = false,),
              Lines    = (linewidth = 2,))

# ---- the figure next to the script that draws it, examples/lorenz/<name>.png ----
finish(fig, name) = save(joinpath(@__DIR__, "$(name).png"), fig; px_per_unit=2)
