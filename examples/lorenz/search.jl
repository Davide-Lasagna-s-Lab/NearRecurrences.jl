# The plain search on the Lorenz system, for the periods from 0.5 to 4 below the distance 0.1, along a
# trajectory of 10⁴ time units. Two runs differ by apart alone, 0.5 and 5 time units, a half and some
# four Lyapunov times (the leading exponent is near 0.9). (a) The matrix of distances of the first
# 60 time units, the minima below the threshold, and those that each rule of the closest keeps.
# (b) The closest recurrence of the wide rule on the attractor, with the seeds of four segments.
#
#     julia --project=examples examples/lorenz/search.jl
include(joinpath(@__DIR__, "lorenz.jl"))


# ---- the three searches on one stream: every minimum, apart 0.5 and apart 5. Returns their
#      recurrences, the closest of the wide rule with its seeds and its stretch, and the samples of
#      the first 60 time units for the matrix ----
function stream(x)
    settings = (spacing   = 0.05,
                windows   = [(0.5, 4.0)],
                threshold = 0.1)
    minima   = NearRecurrences.Recurrences(x, distance; settings..., apart = 0.0)
    narrow   = NearRecurrences.Recurrences(x, distance; settings..., apart = 0.5)
    wide     = NearRecurrences.Recurrences(x, distance; settings..., apart = 5.0)

    found   = ([], [], [])
    closest = nothing
    first60 = Vector{Float64}[]
    for k in 1:200_000
        advance!(x)
        k <= 1200 &&
            push!(first60, copy(x))

        append!(found[1], push!(minima, x))
        append!(found[2], push!(narrow, x))

        # ---- the wide rule, and the closest with its seeds and its stretch, while they are held ----
        for r in push!(wide, x)
            push!(found[3], r)
            if closest === nothing || r.d < closest.r.d
                closest = (r       = r,
                           seeds   = NearRecurrences.seeds(wide, r, 4),
                           stretch = NearRecurrences.seeds(wide, r, r.lag).states)
            end
        end
    end
    for (R, list) in zip((minima, narrow, wide), found)
        append!(list, NearRecurrences.finish!(R))
    end

    return found, closest, first60, minima.lags[1]
end

found, closest, first60, lags = stream(start())

@printf("minima %d, kept with apart 0.5: %d, with apart 5: %d, in 10⁴ time units\n", length.(found)...)
@printf("the closest: t = %.2f, T = %.2f, d = %.1e, seeds at the samples %s\n",
        closest.r.t, closest.r.T, closest.r.d, string(closest.seeds.offsets))


# ---- the figure ----
fig = Figure(size=(1300, 470))

# ---- (a) log10 d over the start and the period, and the three sets of minima ----
n    = length(first60) - last(lags)
D    = [log10(distance(first60[i], first60[i + ℓ])) for i in 1:n, ℓ in lags]
seen(list) = [r for r in list if r.start <= n]

ax = Axis(fig[1, 1]; title="(a) the matrix of distances, the minima and the rule of the closest",
          xlabel="start t", ylabel="period T")
hm = heatmap!(ax, 0.05 .* (0:(n - 1)), 0.05 .* lags, D;
              colormap=cgrad(["#0D47A1", "#1976D2", "#90CAF9", "#FFFFFF"]), colorrange=(-2, 0))
scatter!(ax, [r.t for r in seen(found[1])], [r.T for r in seen(found[1])]; color=:white,
         strokecolor=GREY, strokewidth=1.5, markersize=9, label="minima below 0.1")
scatter!(ax, [r.t for r in seen(found[2])], [r.T for r in seen(found[2])]; color=PURPLE,
         markersize=11, marker=:rect, label="kept, apart 0.5")
scatter!(ax, [r.t for r in seen(found[3])], [r.T for r in seen(found[3])]; color=RED,
         strokecolor=:white, strokewidth=1.5, markersize=17, marker=:star5, label="kept, apart 5")
Colorbar(fig[1, 2], hm; label="log10 d", ticks=(-2:0, ["-2", "-1", "0"]))
Legend(fig[2, 1], ax; orientation=:horizontal)

# ---- (b) the closest on the attractor, the plane (x, z) ----
s  = fine(closest.stretch[1], 5*closest.r.lag)
A  = fine(first60[1], 5*length(first60))
g  = closest.seeds.states
ax = Axis(fig[1, 3]; title=@sprintf("(b) the closest, T = %.2f", closest.r.T), xlabel="x",
          ylabel="z", aspect=1)
lines!(ax, [u[1] for u in A], [u[3] for u in A]; color=(GREY, 0.25), linewidth=1)
lines!(ax, [u[1] for u in s], [u[3] for u in s]; color=:black, linewidth=2)
scatter!(ax, [u[1] for u in g[2:4]], [u[3] for u in g[2:4]]; color=PURPLE, markersize=15,
         label="seeds")
scatter!(ax, [g[5][1]], [g[5][3]]; color=RED, markersize=11, label="end")
scatter!(ax, [g[1][1]], [g[1][3]]; color=:transparent, strokecolor=BLUE, strokewidth=3, markersize=24,
         label="start")
Legend(fig[2, 3], ax; orientation=:horizontal)

colsize!(fig.layout, 1, Relative(0.6))
finish(fig, "search")
