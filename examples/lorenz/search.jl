# The plain search on the Lorenz system, for the periods from 0.5 to 4 below the distance 0.1, along a
# trajectory of 10⁴ time units, with apart_t of 5 time units, some four Lyapunov times (the leading
# exponent is near 0.9), and the whole window along the period. (a) The matrix of distances of the
# first 60 time units, with the recurrences. (b) The closest recurrence on the attractor, with the
# seeds of four segments.
#
#     julia --project=examples examples/lorenz/search.jl
include(joinpath(@__DIR__, "lorenz.jl"))


# ---- the search on one stream. Returns its recurrences, the closest with its seeds and its
#      stretch, and the samples of the first 60 time units for the matrix ----
function stream(x)
    R = NearRecurrences.Recurrences(x, distance;
                                    spacing   = 0.05,
                                    windows   = [(0.5, 4.0)],
                                    threshold = 0.1,
                                    apart_t   = 5.0)

    found   = []
    closest = nothing
    first60 = Vector{Float64}[]
    for k in 1:200_000
        advance!(x)
        k <= 1200 &&
            push!(first60, copy(x))

        # ---- the recurrences, and the closest with its seeds and its stretch, while they are held ----
        for r in push!(R, x)
            push!(found, r)
            if closest === nothing || r.d < closest.r.d
                closest = (r       = r,
                           seeds   = NearRecurrences.seeds(R, r, 4),
                           stretch = NearRecurrences.seeds(R, r, r.lag).states)
            end
        end
    end
    append!(found, NearRecurrences.finish!(R))

    return found, closest, first60, R.lags[1]
end

found, closest, first60, lags = stream(start())

@printf("%d recurrences in 10⁴ time units\n", length(found))
@printf("the closest: t = %.2f, T = %.2f, d = %.1e, seeds at the samples %s\n",
        closest.r.t, closest.r.T, closest.r.d, string(closest.seeds.offsets))


# ---- the figure ----
fig = Figure(size=(1300, 470))

# ---- (a) log10 d over the start and the period, and the recurrences ----
n    = length(first60) - last(lags)
D    = [log10(distance(first60[i], first60[i + ℓ])) for i in 1:n, ℓ in lags]
seen = [r for r in found if r.start <= n]

ax = Axis(fig[1, 1]; title="(a) the matrix of distances and the recurrences",
          xlabel="start t", ylabel="period T")
hm = heatmap!(ax, 0.05 .* (0:(n - 1)), 0.05 .* lags, D;
              colormap=cgrad(["#0D47A1", "#1976D2", "#90CAF9", "#FFFFFF"]), colorrange=(-2, 0))
scatter!(ax, [r.t for r in seen], [r.T for r in seen]; color=RED, strokecolor=:white,
         strokewidth=1.5, markersize=17, marker=:star5, label="recurrences")
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
