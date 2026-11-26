# The probability of a distance on the Lorenz system, with the rotation. The histograms of the search
# give the fraction C(d) of the pairs closer than d, at the lags of the window, exactly at the edges
# of their bins. (a) C against d. (b) The distance against the period of the recurrences, with the
# distances d_p where C reaches the probabilities p, as the note draws them for Kolmogorov flow.
# Periods from 0.5 to 4 below 0.1, apart 5, along 10⁵ time units.
#
#     julia --project=examples examples/lorenz/probability.jl
include(joinpath(@__DIR__, "lorenz.jl"))


# ---- the search on one stream. Returns the recurrences and the search, for its histograms ----
function stream(x)
    R = NearRecurrences.Recurrences(x, rotated;
                                    spacing   = 0.05,
                                    windows   = [(0.5, 4.0)],
                                    threshold = 0.1,
                                    apart     = 5.0,
                                    invariant = invariant)
    found = []
    for _ in 1:2_000_000
        advance!(x)
        append!(found, push!(R, x))
    end
    append!(found, NearRecurrences.finish!(R))

    return found, R
end

found, R = stream(start())

# ---- C at the edges of the bins, and the distances of the probabilities 10^-2, 10^-2.5 and 10^-3,
#      the first edges where C reaches them, all above the first edge, where C is resolved ----
edges = R.edges[1][2:end]
C     = [NearRecurrences.probability(R, exp10(e)) for e in edges]
p     = [-2.0, -2.5, -3.0]
dp    = [exp10(edges[findfirst(>=(exp10(k)), C)]) for k in p]
@assert all(C[1] .< exp10.(p))
@printf("%d recurrences, %d pairs, C(0.1) = %.2e\n", length(found), sum(R.pairs), C[end])
for (k, d) in zip(p, dp)
    @printf("p = 10^%.1f at d = %.4f\n", k, d)
end


# ---- the figure ----
fig    = Figure(size=(1300, 520))
colour = (BLUE, RED, PURPLE)

# ---- (a) log10 C against log10 d ----
keep = C .> 0
ax   = Axis(fig[1, 1]; title="(a) the probability of a distance", xlabel="log10 d",
            ylabel="log10 C(d)", xticks=-2:-1, yticks=-9:0)
lines!(ax, edges[keep], log10.(C[keep]); color=:black, linewidth=2.5)
for (j, k) in enumerate(p)
    scatter!(ax, [log10(dp[j])], [log10(C[findfirst(>=(exp10(k)), C)])]; color=colour[j], markersize=15)
end

# ---- (b) the distance against the period, and the levels ----
ax = Axis(fig[1, 2]; title="(b) the recurrences and the distances of three probabilities",
          xlabel="period T", ylabel="distance d")
scatter!(ax, [r.T for r in found], [r.d for r in found]; color=(GREY, 0.5), markersize=5)
for (j, k) in enumerate(p)
    hlines!(ax, [dp[j]]; color=colour[j], linestyle=:dash, linewidth=2.5, label="p = $(round(exp10(k); sigdigits=1))")
end
Legend(fig[1, 3], ax)

colsize!(fig.layout, 1, Relative(0.38))
finish(fig, "probability")
