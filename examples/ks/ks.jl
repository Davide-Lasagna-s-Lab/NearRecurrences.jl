# The near recurrences of the Kuramoto–Sivashinsky equation u_t = -u u_x - u_xx - u_xxxx on a periodic
# domain of length 22, with its two symmetries, the shift along x, continuous, and the reflection
# u(x) → -u(-x), discrete. The system is that of the tests, test/systems.jl: the modes k = 1 to 15 of u,
# the products on 48 points, the classical Runge–Kutta scheme on e^{-Λt} a with the linear part Λ taken
# exactly, steps of 0.025, a sample every 0.25. The distance is minimised over the shift and the
# reflection and returns (d, σ, r), and the invariant is the weighted moduli √2 |a_k|.
#
# A stream of 5 × 10⁴ time units, for the periods from 8 to 40 below 0.15, with apart_t 90, about four
# Lyapunov times at the leading exponent near 0.04. (search.png) The distance against the period of
# the recurrences, with the distances where the probability C(d) reaches three levels, and the
# space–time picture of the stretch of the closest recurrence, with the seeds of four segments.
# (bound.png) The bound against the distance for the pairs of every tenth start of 2000 time units,
# with the fraction of the pairs that the stream skips.
#
#     julia --project=examples examples/ks/ks.jl
import NearRecurrences
using CairoMakie
using Printf

include(joinpath(@__DIR__, "..", "..", "test", "systems.jl"))


# ---------------------------------------------------------------------------- #
# the stream                                                                   #
# ---------------------------------------------------------------------------- #

# The stream, the calls of the distance counted. Returns the recurrences, the search, the calls, the
# closest recurrence with its seeds and its stretch, and the samples of the first 2000 time units.
function stream()
    u    = zeros(30)
    u[1] = 0.3
    u[4] = -0.2
    u[7] = 0.1
    transient(ks_advance!, u, 1000)

    calls = Ref(0)
    R     = NearRecurrences.Recurrences(u, (a, b) -> (calls[] += 1; ks_distance(a, b));
                                        spacing   = 0.25,
                                        windows   = [(8.0, 40.0)],
                                        threshold = 0.15,
                                        apart_t   = 90.0,
                                        invariant = ks_invariant)

    found   = []
    closest = nothing
    first2k = Vector{Float64}[]
    for k in 1:200_000
        ks_advance!(u)
        k <= 8000 &&
            push!(first2k, copy(u))
        for r in push!(R, u)
            push!(found, r)
            if closest === nothing || r.d < closest.r.d
                closest = (r       = r,
                           seeds   = NearRecurrences.seeds(R, r, 4),
                           stretch = NearRecurrences.seeds(R, r, r.lag).states)
            end
        end
    end
    append!(found, NearRecurrences.finish!(R))

    return found, R, calls[], closest, first2k
end

found, R, calls, closest, first2k = stream()
pairs = sum(R.pairs)
@printf("%d recurrences in 5 × 10⁴ time units, %d with the reflection, %d pairs, %.1f per cent compared in full\n",
        length(found), count(r -> r.shifts[2] == 1, found), pairs, 100calls/pairs)
@printf("the closest: T = %.2f, d = %.1e, σ = %.3f, r = %d, seeds at the samples %s\n",
        closest.r.T, closest.r.d, closest.r.shifts..., string(closest.seeds.offsets))

# ---- C at the edges of the bins up to the threshold, and the distances of three levels of it, half
#      a decade apart ----
edges  = [e for e in R.edges[1][2:end] if exp10(e) <= 0.15]
C      = [NearRecurrences.probability(R, exp10(e)) for e in edges]
top    = floor(2*log10(C[end]))/2
levels = [k for k in top:-0.5:(top - 1.0) if C[1] < exp10(k)]
dp     = [exp10(edges[findfirst(>=(exp10(k)), C)]) for k in levels]
for (k, d) in zip(levels, dp)
    @printf("C = 10^%.1f at d = %.4f\n", k, d)
end

# ---- the distance and the bound of the pairs of every tenth start of the first 2000 time units ----
function sampled(samples)
    d = Float64[]
    β = Float64[]
    for i in 1:10:(length(samples) - 160), ℓ in 32:160
        a = samples[i]
        b = samples[i + ℓ]
        push!(d, first(ks_distance(a, b)))
        push!(β, sqrt(sum(abs2, ks_invariant(b) .- ks_invariant(a)))/sqrt(sum(abs2, ks_invariant(a))))
    end
    return d, β
end

d, β = sampled(first2k)
@printf("%d pairs, largest β/d %.4f\n", length(d), maximum(β ./ d))


# ---------------------------------------------------------------------------- #
# the figures                                                                  #
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
              Legend   = (framevisible = false,))

# ---- u on 64 points of the domain, from the modes of the state ----
const X = range(0, KS_L; length=65)[1:64]
field(u) = [sum(2*real(complex(u[2k - 1], u[2k])*cis(KS_Q[k]*x)) for k in 1:KS_K) for x in X]

# ---- search.png: (a) the recurrences and the levels, (b) the stretch of the closest ----
fig    = Figure(size=(1300, 560))
colour = (BLUE, RED, PURPLE)

ax = Axis(fig[1, 1]; title="(a) the recurrences, by the reflection that closes them",
          xlabel="period T", ylabel="distance d")
none = [r for r in found if r.shifts[2] == 0]
once = [r for r in found if r.shifts[2] == 1]
scatter!(ax, [r.T for r in none], [r.d for r in none]; color=(GREY, 0.6), markersize=8,
         label="r = 0, the end shifted, $(length(none))")
scatter!(ax, [r.T for r in once], [r.d for r in once]; color=(:black, 0.8), markersize=10,
         marker=:utriangle, label="r = 1, and reflected, $(length(once))")
for (j, k) in enumerate(levels)
    hlines!(ax, [dp[j]]; color=colour[j], linestyle=:dash, linewidth=2.5, label="p = $(round(exp10(k); sigdigits=1))")
end
Legend(fig[2, 1], ax; orientation=:horizontal, nbanks=2)

s  = closest.stretch
t  = 0.25 .* (0:(length(s) - 1))
U  = reduce(hcat, field.(s))
ax = Axis(fig[1, 2]; title=@sprintf("(b) the closest, T = %.2f, with its seeds", closest.r.T),
          xlabel="x", ylabel="time from the start")
hm = heatmap!(ax, X, t, U; colormap=cgrad(["#0D47A1", "#1976D2", "#FFFFFF", "#E53935", "#B71C1C"]),
              colorrange=(-3, 3))
hlines!(ax, 0.25 .* closest.seeds.offsets; color=PURPLE, linewidth=2.5)
Colorbar(fig[1, 3], hm; label="u")

colsize!(fig.layout, 1, Relative(0.55))
save(joinpath(@__DIR__, "search.png"), fig; px_per_unit=2)

# ---- bound.png: log10 β against log10 d, as a density of the pairs ----
fig     = Figure(size=(760, 700))
skipped = 1 - calls/pairs
bins    = -1.6:0.02:0.4
counts  = zeros(length(bins) - 1, length(bins) - 1)
for (x, y) in zip(log10.(d), log10.(β))
    i = searchsortedlast(bins, x)
    j = searchsortedlast(bins, y)
    1 <= i < length(bins) && 1 <= j < length(bins) &&
        (counts[i, j] += 1)
end
centres = bins[1:(end - 1)] .+ 0.01

ax = Axis(fig[1, 1]; title="the bound against the distance, for the pairs of a stream",
          xlabel="log10 d, the distance up to the shift and the reflection", ylabel="log10 β, the bound",
          xticks=-1:0, yticks=-1:0, limits=((-1.6, 0.4), (-1.6, 0.4)), aspect=1)
poly!(ax, Point2f[(-1.6, log10(0.15)), (0.4, log10(0.15)), (0.4, 0.4), (-1.6, 0.4)]; color=(GREY, 0.12))
heatmap!(ax, centres, centres, map(c -> c > 0 ? log10(c) : NaN, counts);
         colormap=cgrad(["#BBDEFB", "#1976D2", "#0D47A1"]))
lines!(ax, [-1.6, 0.4], [-1.6, 0.4]; color=:black, linewidth=1.5, label="β = d")
hlines!(ax, [log10(0.15)]; color=RED, linestyle=:dash, linewidth=2, label="the threshold 0.15")
vlines!(ax, [log10(0.15)]; color=RED, linestyle=:dash, linewidth=2)
text!(ax, -1.55, 0.35; text=@sprintf("skipped, β ≥ 0.15,\n%.0f %% of the pairs", 100skipped),
      align=(:left, :top), color=GREY, fontsize=17)
axislegend(ax; position=:rb)

save(joinpath(@__DIR__, "bound.png"), fig; px_per_unit=2)
println("written examples/ks/search.png and bound.png")
