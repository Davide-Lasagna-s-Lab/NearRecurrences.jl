# The search with the symmetry of the Lorenz system, the rotation R(x, y, z) = (-x, -y, z) by π about
# the z axis. The distance is the smaller of ‖b - a‖ and ‖R b - a‖ over ‖a‖ and returns the rotations
# r of b that reach it, so that a recurrence may close after the rotation, as the two halves of a
# symmetric orbit do. This is the small model of the shift and the translation of Kolmogorov flow.
# Periods from 0.5 to 4 below 0.1, apart 5, along 10⁴ time units. (a) The distance against the
# period of the recurrences, by their rotation. (b) The closest recurrence that closes after the
# rotation, its end rotated back near its start.
#
#     julia --project=examples examples/lorenz/symmetry.jl
include(joinpath(@__DIR__, "lorenz.jl"))


# ---- the search on one stream. Returns the recurrences, the closest of those with r = 1 with its
#      stretch, and the samples of the first 60 time units for the attractor ----
function stream(x)
    R = NearRecurrences.Recurrences(x, rotated;
                                    spacing   = 0.05,
                                    windows   = [(0.5, 4.0)],
                                    threshold = 0.1,
                                    apart     = 5.0)
    found   = []
    shown   = nothing
    first60 = Vector{Float64}[]
    for k in 1:200_000
        advance!(x)
        k <= 1200 &&
            push!(first60, copy(x))

        for r in push!(R, x)
            push!(found, r)
            if r.shifts[1] == 1 && (shown === nothing || r.d < shown.r.d)
                shown = (r       = r,
                         stretch = NearRecurrences.seeds(R, r, r.lag).states)
            end
        end
    end
    append!(found, NearRecurrences.finish!(R))

    return found, shown, first60
end

found, shown, first60 = stream(start())
once = [r for r in found if r.shifts[1] == 1]
none = [r for r in found if r.shifts[1] == 0]
@printf("%d recurrences, %d of them closing after the rotation\n", length(found), length(once))
@printf("the closest of these: T = %.2f, d = %.1e\n", shown.r.T, shown.r.d)


# ---- the figure ----
fig = Figure(size=(1300, 520))

# ---- (a) the distance against the period, by the rotation that closes the recurrence ----
ax = Axis(fig[1, 1]; title="(a) the recurrences, by the rotation that closes them",
          xlabel="period T", ylabel="distance d")
scatter!(ax, [r.T for r in none], [r.d for r in none]; color=(BLUE, 0.7), markersize=9,
         label="r = 0, the end itself, $(length(none))")
scatter!(ax, [r.T for r in once], [r.d for r in once]; color=(RED, 0.7), markersize=11,
         marker=:utriangle, label="r = 1, the end rotated, $(length(once))")
Legend(fig[2, 1], ax; orientation=:horizontal)

# ---- (b) the closest that closes after the rotation, in the plane (x, z), where R is the mirror
#      x → -x ----
s  = fine(shown.stretch[1], 5*shown.r.lag)
A  = fine(first60[1], 5*length(first60))
a  = s[1]
b  = s[end]
ax = Axis(fig[1, 2]; title=@sprintf("(b) one that closes after the rotation, T = %.2f", shown.r.T),
          xlabel="x", ylabel="z", aspect=1)
lines!(ax, [u[1] for u in A], [u[3] for u in A]; color=(GREY, 0.25), linewidth=1)
vlines!(ax, [0.0]; color=GREY, linestyle=:dash, linewidth=1.5)
lines!(ax, [u[1] for u in s], [u[3] for u in s]; color=:black, linewidth=2)
lines!(ax, [b[1], -b[1]], [b[3], b[3]]; color=PURPLE, linestyle=:dot, linewidth=2)
scatter!(ax, [b[1]], [b[3]]; color=RED, markersize=15, label="end")
scatter!(ax, [-b[1]], [b[3]]; color=PURPLE, markersize=13, marker=:diamond, label="R of the end")
scatter!(ax, [a[1]], [a[3]]; color=:transparent, strokecolor=BLUE, strokewidth=3, markersize=26,
         label="start")
Legend(fig[2, 2], ax; orientation=:horizontal)

colsize!(fig.layout, 1, Relative(0.5))
finish(fig, "symmetry")
