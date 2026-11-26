# The figure of bound.jl, from bound.csv: (a) the bound against the distance of the pairs of Kolmogorov
# flow as a density, with the threshold 0.118 and the fraction of the pairs that the stream skips, and
# (b) the fraction of the pairs compared in full against the threshold.
#
#     julia --project=examples examples/kolmogorov/figure.jl
using CairoMakie
using Printf


# ---- the pairs, and the counts of the stream from the first line ----
table   = readlines(joinpath(@__DIR__, "bound.csv"))
counted = parse.(Int, match(r"pairs (\d+), compared (\d+)", table[1]).captures)
pairs   = counted[1]
calls   = counted[2]
rows    = [parse.(Float64, split(l, ",")) for l in table[3:end]]
d       = [r[1] for r in rows]
β       = [r[2] for r in rows]


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
              Legend   = (framevisible = false,))

fig       = Figure(size=(1300, 600))
threshold = 0.118
skipped   = 1 - calls/pairs

# ---- (a) the density of the pairs over log10 d and log10 β, bins of 0.02 ----
edges  = -1.6:0.02:0.4
counts = zeros(length(edges) - 1, length(edges) - 1)
for (x, y) in zip(log10.(d), log10.(β))
    i = searchsortedlast(edges, x)
    j = searchsortedlast(edges, y)
    1 <= i < length(edges) && 1 <= j < length(edges) &&
        (counts[i, j] += 1)
end
centres = edges[1:(end - 1)] .+ 0.01

ax = Axis(fig[1, 1]; title="(a) the bound against the distance, for the pairs of a stream",
          xlabel="log10 d, the distance up to the shift and the translations",
          ylabel="log10 β, the bound", xticks=-1:0, yticks=-1:0,
          limits=((-1.6, 0.4), (-1.6, 0.4)), aspect=1)
poly!(ax, Point2f[(-1.6, log10(threshold)), (0.4, log10(threshold)), (0.4, 0.4), (-1.6, 0.4)];
      color=(GREY, 0.12))
hm = heatmap!(ax, centres, centres, map(c -> c > 0 ? log10(c) : NaN, counts);
              colormap=cgrad(["#BBDEFB", "#1976D2", "#0D47A1"]))
lines!(ax, [-1.6, 0.4], [-1.6, 0.4]; color=:black, linewidth=1.5, label="β = d")
hlines!(ax, [log10(threshold)]; color=RED, linestyle=:dash, linewidth=2, label="the threshold 0.118")
vlines!(ax, [log10(threshold)]; color=RED, linestyle=:dash, linewidth=2)
text!(ax, -1.55, 0.35; text=@sprintf("skipped, β ≥ 0.118,\n%.3f %% of the pairs", 100skipped),
      align=(:left, :top), color=GREY, fontsize=17)
axislegend(ax; position=:rb)
Colorbar(fig[1, 2], hm; label="log10 pairs per bin")

# ---- (b) the fraction of the pairs compared in full against the threshold ----
thresholds = range(0.05, 0.4; length=36)
compared   = [count(<(t), β)/length(β) for t in thresholds]
keep       = compared .> 0
ax = Axis(fig[1, 3]; title="(b) the pairs compared in full", xlabel="threshold",
          ylabel="log10 fraction of the pairs", yticks=-5:0)
lines!(ax, thresholds[keep], log10.(compared[keep]); color=:black, linewidth=2.5)
vlines!(ax, [threshold]; color=RED, linestyle=:dash, linewidth=2)
scatter!(ax, [threshold], [log10(calls/pairs)]; color=RED, markersize=15)

colsize!(fig.layout, 1, Aspect(1, 1.0))
colsize!(fig.layout, 3, Aspect(1, 0.8))
resize_to_layout!(fig)
save(joinpath(@__DIR__, "bound.png"), fig; px_per_unit=2)
println("written examples/kolmogorov/bound.png")
