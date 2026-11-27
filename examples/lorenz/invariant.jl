# The invariant bound on the search with the rotation of the Lorenz system. The invariant
# v = (√(x² + y²), z) does not change under the rotation, has the norm of the state, and does not
# stretch distances, so that β = ‖v(b) - v(a)‖/‖v(a)‖ is at most the distance d minimised over the
# rotation. A pair whose β is not below the threshold is not compared in full. Periods from 0.5 to 4
# below 0.1, apart_t 5, along 10⁴ time units, with the invariant and without it, the calls of the
# distance counted. The figure draws β against d for the pairs of every tenth start of the first
# 1000 time units, all below the diagonal, with the threshold and the pairs it skips.
#
#     julia --project=examples examples/lorenz/invariant.jl
include(joinpath(@__DIR__, "lorenz.jl"))


# ---- the two searches on one stream, with the calls of the distance counted. Returns their
#      recurrences, the calls, the pairs compared, and the samples of the first 1000 time units ----
function stream(x)
    calls    = zeros(Int, 2)
    settings = (spacing   = 0.05,
                windows   = [(0.5, 4.0)],
                threshold = 0.1,
                apart_t   = 5.0)
    full     = NearRecurrences.Recurrences(x, (a, b) -> (calls[1] += 1; rotated(a, b)); settings...)
    bound    = NearRecurrences.Recurrences(x, (a, b) -> (calls[2] += 1; rotated(a, b)); settings...,
                                           invariant = invariant)
    found    = ([], [])
    first1000 = Vector{Float64}[]
    for k in 1:200_000
        advance!(x)
        k <= 20_000 &&
            push!(first1000, copy(x))
        append!(found[1], push!(full, x))
        append!(found[2], push!(bound, x))
    end
    append!(found[1], NearRecurrences.finish!(full))
    append!(found[2], NearRecurrences.finish!(bound))

    return found, calls, sum(full.pairs), first1000, full.lags[1]
end

found, calls, pairs, first1000, lags = stream(start())
skipped = 1 - calls[2]/pairs
@printf("the same recurrences with the invariant: %s, %d of them\n", found[1] == found[2], length(found[1]))
@printf("calls of the distance: %d without the invariant, %d with it, for %d pairs: %.1f per cent skipped\n",
        calls[1], calls[2], pairs, 100skipped)


# ---- the bound and the distance of the pairs of every tenth start ----
d = Float64[]
β = Float64[]
for i in 1:10:(length(first1000) - last(lags)), ℓ in lags
    a = first1000[i]
    b = first1000[i + ℓ]
    push!(d, first(rotated(a, b)))
    push!(β, sqrt(sum(abs2, invariant(b) .- invariant(a)))/sqrt(sum(abs2, invariant(a))))
end
@printf("largest β/d over %d pairs: %.4f\n", length(d), maximum(β ./ d))


# ---- the figure: log10 β against log10 d ----
fig = Figure(size=(760, 700))
ax  = Axis(fig[1, 1]; title="the bound against the distance, for the pairs of a stream",
           xlabel="log10 d, the distance up to the rotation", ylabel="log10 β, the bound",
           xticks=-3:0, yticks=-4:0, limits=((-3, 0.3), (-4, 0.3)), aspect=1)

# the pairs skipped, β at or above the threshold
poly!(ax, Point2f[(-1, -1), (0.3, -1), (0.3, 0.3), (-1, 0.3)]; color=(GREY, 0.15))
text!(ax, -0.95, 0.2; text=@sprintf("skipped, β ≥ 0.1,\n%.0f %% of the pairs", 100skipped),
      align=(:left, :top), color=GREY, fontsize=17)

scatter!(ax, log10.(d), log10.(β); color=(BLUE, 0.25), markersize=3, strokewidth=0)
lines!(ax, [-3, 0.3], [-3, 0.3]; color=:black, linewidth=1.5, label="β = d")
hlines!(ax, [-1]; color=RED, linestyle=:dash, linewidth=2, label="the threshold 0.1")
vlines!(ax, [-1]; color=RED, linestyle=:dash, linewidth=2)
axislegend(ax; position=:lt)

finish(fig, "invariant")
