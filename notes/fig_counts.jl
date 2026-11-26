# The figure of the counts of the note (recurrence-search.tex, figure fig:numbers) and the numbers of
# its last section, from the tables of data/.
#
# (a) The distance d of every recurrence of Kolmogorov flow found with the shift along x alone against
#     its period T, with the distances of the probabilities 1e-5, 1e-6 and 1e-7.
# (b) The recurrences per 1e5 time units of simulation in a band of periods T ± 5 % below the same
#     distances, with the shift along x alone, filled, and with the translations along y as well,
#     hollow, whose windows of periods are scaled to a band of ± 5 % about their middle. The bars are
#     ± one standard deviation of a Poisson count, √n.
#
# The tables:
#     data/shift-x.csv         the recurrences with the shift along x alone, by band
#     data/shift-x-runs.csv    their runs, the bands searched, the time simulated and the costs
#     data/shift-xy.csv        the recurrences with the translations along y, by window
#     data/shift-xy-runs.csv   their runs, windows, thresholds and times
#     data/bound-runs.csv      the runs with the bound by the moduli, their costs and the fraction of
#                              the pairs compared in full
#
# Writes figures/counts.pdf and .png, and prints the counts, the exponent μ, the gain from the
# translations and the costs.
#
#     julia --project=<an environment with CairoMakie> notes/fig_counts.jl
using CairoMakie
using Printf
using Statistics


# ---- places and settings ----
const OUT   = joinpath(@__DIR__, "figures", "counts")
const DATA  = joinpath(@__DIR__, "data")
const D     = (0.1592, 0.118, 0.091)                    # the distances of p = 1e-5, 1e-6 and 1e-7
const BANDS = [20, 30, 50, 70, 100, 150, 200, 300, 500, 1000, 2000]

# ---- the look: blue, red, purple and grey, the fonts of LaTeX, no top and right spines ----
const BLUE   = Makie.to_color("#1976D2")
const RED    = Makie.to_color("#E53935")
const PURPLE = Makie.to_color("#5E35B1")
const GREY   = Makie.to_color("#757575")

set_theme!(theme_latexfonts())
update_theme!(fontsize       = 12,
              figure_padding = (2, 6, 2, 2),
              Axis           = (topspinevisible   = false,
                                rightspinevisible = false,
                                xgridcolor        = RGBAf(0, 0, 0, 0.08),
                                ygridcolor        = RGBAf(0, 0, 0, 0.08),
                                xticklabelsize    = 10,
                                yticklabelsize    = 10))


# ---- the tables, as vectors of rows split at the commas ----
table(name) = [split(l, ",") for l in readlines(joinpath(DATA, name))[2:end] if !isempty(l)]

alone  = table("shift-x.csv")
runs   = table("shift-x-runs.csv")
both   = table("shift-xy.csv")
bruns  = table("shift-xy-runs.csv")
bound  = table("bound-runs.csv")


# ---- the shift along x alone: the time simulated for each band, the counts below each distance ----
band  = [parse(Int, r[1]) for r in alone]
dist  = [parse(Float64, r[5]) for r in alone]
simul = Dict(b => 0.0 for b in BANDS)
for r in runs, b in parse.(Int, split(r[2], "/"))
    simul[b] += parse(Float64, r[3])
end
count(b, d) = sum((band .== b) .& (dist .< d))
rate(b, d)  = 1e5*count(b, d)/simul[b]

println("time units ", sum(parse(Float64, r[3]) for r in runs), ", recurrences ", length(dist))
println("band, time units, counts below p = 1e-5, 1e-6, 1e-7, rates per 1e5 time units")
for b in BANDS
    @printf "%5d %9.3g  %5d %4d %3d  %7.3f %7.3f %7.4f\n" b simul[b] count(b, D[1]) count(b, D[2]) count(b, D[3]) rate(b, D[1]) rate(b, D[2]) rate(b, D[3])
end

# the exponent μ of the counts below d, N(< d) ∝ d^μ, by maximum likelihood on those below D[1]
μ(d)  = length(d)/sum(log.(D[1] ./ d))
below = dist[dist .< D[1]]
@printf "μ pooled %.2f, in the bands 20 and 30 %.2f; fractions below 1e-6 and 1e-7: %.3f %.3f\n" μ(below) μ(dist[(band .<= 30) .& (dist .< D[1])]) mean(below .< D[2]) mean(below .< D[3])


# ---- with the translations along y: rates per 1e5 time units by window, and with no translation ----
# a window is the pair (study, lo:hi), since one window ran in two studies with two thresholds
window = Dict{Tuple{String, String}, Vector{Float64}}()   # the distances by window
level  = Dict{Tuple{String, String}, Vector{Float64}}()   # those with no translation
for r in both
    w = (String(r[1]), String(r[2]))
    push!(get!(window, w, Float64[]), parse(Float64, r[6]))
    parse(Int, r[5]) == 0 &&
        push!(get!(level, w, Float64[]), parse(Float64, r[6]))
end
otime = Dict{Tuple{String, String}, Float64}()             # the time simulated by window
ocut  = Dict{Tuple{String, String}, Float64}()             # the threshold of the window
for r in bruns, lohi in split(r[3], " ")
    w        = (String(r[1]), String(lohi))
    otime[w] = get(otime, w, 0.0) + parse(Float64, r[6])
    ocut[w]  = parse(Float64, r[4])
end
centre(w) = mean(parse.(Float64, split(w[2], ":")))
width(w)  = -(reverse(parse.(Float64, split(w[2], ":")))...)
orate(w, d; v=window) = 1e5*sum(get(v, w, Float64[]) .< d)/otime[w]*(0.1centre(w)/width(w))

gain = sum(sum(window[w] .< 0.15) for w in keys(window))/sum(sum(level[w] .< 0.15) for w in keys(level))
@printf "gain from the translations along y, all windows below 0.15: %.2f\n" gain


# ---- the costs per time unit ----
println("shift along x: bands, CPU, runs, median s per time unit of integration and per distance, pairs per sample")
for k in unique((r[2], r[8]) for r in runs)
    rs = [r for r in runs if (r[2], r[8]) == k]
    @printf "%-28s %-10s %3d  %.4f %.2e  %s\n" k[1] k[2] length(rs) median(parse(Float64, r[5]) for r in rs) median(parse(Float64, r[6]) for r in rs) rs[1][7]
end
println("with the bound, EPYC 7502: windows, runs, median s per time unit of integration and of distances, s per pair, fraction in full")
for w in unique(r[2] for r in bound)
    rs = [r for r in bound if r[2] == w]
    @printf "%-12s %2d  %.4f %.4f  %.2e  %.1e\n" first(w, 12) length(rs) median(parse(Float64, r[4]) for r in rs) median(parse(Float64, r[5]) for r in rs) median(parse(Float64, r[5])/(10parse(Float64, r[6])) for r in rs) median(parse(Float64, r[7]) for r in rs)
end


# ---- the figure ----
fig    = Figure(size=(624, 205))
bticks = ([20, 50, 100, 200, 500, 1000, 2000, 5000, 10000], ["20", "50", "100", "200", "500", "1000", "2000", "5000", "10⁴"])
colour = (BLUE, RED, PURPLE)
marker = (:circle, :rect, :utriangle)
labels = (L"p = 10^{-5}", L"p = 10^{-6}", L"p = 10^{-7}")

# ---- (a) the distance against the period, coloured by the probability below which it falls ----
ax = Axis(fig[1, 1]; xscale=log2, xticks=(bticks[1][1:7], bticks[2][1:7]), xlabel=L"\text{period } T",
          ylabel=L"\text{distance } d", limits=((17, 2300), (0.02, 0.17)))
class = [d < D[3] ? 3 : d < D[2] ? 2 : 1 for d in dist]
T     = [parse(Float64, r[4]) for r in alone]
for k in 1:3
    scatter!(ax, T[class .== k], dist[class .== k]; color=(colour[k], 0.6), markersize=2.5, strokewidth=0)
end
for k in 1:3
    hlines!(ax, [D[k]]; color=colour[k], linewidth=0.9, linestyle=:dash)
    text!(ax, 740, D[k]; text=labels[k], align=(:center, :bottom), fontsize=8, color=colour[k])
end

# ---- (b) log10 of the recurrences per 1e5 time units in a band of ± 5 %, with integer ticks ----
ax = Axis(fig[1, 2]; xscale=log2, xticks=bticks, xlabel=L"\text{period } T", yticks=-2:3,
          ylabel=L"\log_{10}\,\text{count per } 10^5 \text{ time units}", limits=((15, 13000), (-2.3, 3.3)))
for k in 1:3
    x = Float64[]
    y = Float64[]
    e = Float64[]
    for b in BANDS
        n = count(b, D[k])
        n > 0 ||
            continue
        push!(x, b)
        push!(y, rate(b, D[k]))
        push!(e, rate(b, D[k])/sqrt(n))
    end
    rangebars!(ax, x, log10.(max.(y .- e, 5.5e-3)), log10.(y .+ e); color=colour[k], linewidth=0.8)
    scatter!(ax, x, log10.(y); color=colour[k], marker=marker[k], markersize=6, label=labels[k])

    # with the translations along y, where the threshold of the window lies above the distance
    ws = [w for w in keys(otime) if ocut[w] > D[k] && orate(w, D[k]) > 0]
    scatter!(ax, centre.(ws), [log10(orate(w, D[k])) for w in ws]; color=:white, strokecolor=colour[k],
             strokewidth=1, marker=marker[k], markersize=6)
end
axislegend(ax; position=:lt, framevisible=false, labelsize=9, patchsize=(10, 8), rowgap=0)
text!(ax, 11000, -2.1; text="filled: shift along x\nhollow: and along y", align=(:right, :bottom), fontsize=9,
      color=GREY)

Label(fig[1, 1, TopLeft()], "(a)"; fontsize=12, padding=(0, 10, 8, 0), halign=:right, valign=:bottom)
Label(fig[1, 2, TopLeft()], "(b)"; fontsize=12, padding=(0, 10, 8, 0), halign=:right, valign=:bottom)
colgap!(fig.layout, 20)
mkpath(dirname(OUT))
save(OUT*".pdf", fig)
save(OUT*".png", fig; px_per_unit=3)
println("written ", OUT, ".pdf")
