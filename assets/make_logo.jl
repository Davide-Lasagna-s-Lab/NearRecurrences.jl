# Draws assets/logo.svg: a stretch of the Lorenz attractor in the plane (x, z) that leaves its start,
# a blue dot, winds once round each wing and comes back near it, a red dot, with the gap between the
# two in purple, the near recurrence. The stretch is one of the recurrences that the package finds
# along 200 time units, for the periods from 1.5 to 1.6, where the trajectory winds once round each
# wing, below the distance 0.25, the one whose gap is the widest in the plane (x, z) drawn, so that
# it shows at the size of the logo.
#
#     julia --project=. assets/make_logo.jl
import NearRecurrences
using Printf


# ---- the Lorenz system at σ = 10, ρ = 28, β = 8/3, one step of 0.01 of the classical scheme ----
lorenz(x) = [10*(x[2] - x[1]), x[1]*(28 - x[3]) - x[2], x[1]*x[2] - 8/3*x[3]]

function step(x)
    k₁ = lorenz(x)
    k₂ = lorenz(x .+ 0.005 .* k₁)
    k₃ = lorenz(x .+ 0.005 .* k₂)
    k₄ = lorenz(x .+ 0.01 .* k₃)
    return x .+ (0.01/6) .* (k₁ .+ 2 .* k₂ .+ 2 .* k₃ .+ k₄)
end

# ---- one sample ahead, five steps, and the distance ‖b - a‖/‖a‖ ----
function advance!(x)
    for _ in 1:5
        x .= step(x)
    end
    return x
end

distance(a, b) = sqrt(sum(abs2, b .- a))/sqrt(sum(abs2, a))


# ---- the recurrence whose gap is the widest in the plane (x, z): its start and its lag in samples ----
function recurrence()
    x = [1.0, 1.0, 1.0]
    for _ in 1:400
        advance!(x)
    end
    R = NearRecurrences.Recurrences(x, distance;
                                    spacing   = 0.05,
                                    windows   = [(1.5, 1.6)],
                                    threshold = 0.25,
                                    apart_t   = 5.0)
    best = (gap = 0.0, start = Float64[], lag = 0)
    for _ in 1:4000
        advance!(x)
        for r in push!(R, x)
            g   = NearRecurrences.seeds(R, r, 1).states
            gap = hypot(g[2][1] - g[1][1], g[2][3] - g[1][3])
            gap > best.gap &&
                (best = (gap = gap, start = copy(g[1]), lag = r.lag))
        end
    end
    return best.start, best.lag
end

start, lag = recurrence()

# ---- the stretch, every step of 0.01 from the start to the end ----
stretch = [start]
for _ in 1:(5lag)
    push!(stretch, step(stretch[end]))
end


# ---- the plane (x, z) into the box of the drawing, 320 wide and 190 high about its centre ----
xs       = [p[1] for p in stretch]
zs       = [p[3] for p in stretch]
scale    = min(320/(maximum(xs) - minimum(xs)), 190/(maximum(zs) - minimum(zs)))
place(p) = (scale*(p[1] - (maximum(xs) + minimum(xs))/2), -scale*(p[3] - (maximum(zs) + minimum(zs))/2))

path = "M" * join([@sprintf("%.1f %.1f", place(p)...) for p in stretch], " L")
a    = place(stretch[1])
b    = place(stretch[end])


# ---- the svg ----

svg = """
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 960 300" role="img" aria-labelledby="title desc">
  <title id="title">NearRecurrences.jl logo</title>
  <desc id="desc">A stretch of the Lorenz attractor that leaves its start, a blue dot, winds once round each wing and comes back near it, a red dot, with the small gap between the two in purple, beside the NearRecurrences.jl wordmark.</desc>
  <defs>
    <linearGradient id="bg" x1="0" x2="1" y1="0" y2="1">
      <stop offset="0" stop-color="#ffffff"/>
      <stop offset="1" stop-color="#f2f5f8"/>
    </linearGradient>
    <filter id="panel-shadow" x="-12%" y="-25%" width="124%" height="150%">
      <feDropShadow dx="0" dy="9" stdDeviation="11" flood-color="#172033" flood-opacity="0.14"/>
    </filter>
    <filter id="soft" x="-30%" y="-30%" width="160%" height="160%">
      <feDropShadow dx="0" dy="3" stdDeviation="3" flood-color="#172033" flood-opacity="0.10"/>
    </filter>
  </defs>

  <!-- Main panel -->
  <rect x="24" y="24" width="912" height="252" rx="24" fill="url(#bg)" filter="url(#panel-shadow)"/>

  <!-- Left illustration: the stretch and its near return -->
  <g transform="translate(222 150)">
    <path d="$path" fill="none" stroke="#757575" stroke-width="4" stroke-linejoin="round" stroke-linecap="round"/>
    <g filter="url(#soft)">
      $(@sprintf("<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#5E35B1\" stroke-width=\"5\" stroke-linecap=\"round\"/>", a..., b...))
      $(@sprintf("<circle cx=\"%.1f\" cy=\"%.1f\" r=\"6.5\" fill=\"#1976D2\"/>", a...))
      $(@sprintf("<circle cx=\"%.1f\" cy=\"%.1f\" r=\"6.5\" fill=\"#E53935\"/>", b...))
    </g>
  </g>

  <!-- Wordmark and subtitle -->
  <g transform="translate(430 96)">
    <text x="0" y="48" font-family="Inter, Avenir, Helvetica, Arial, sans-serif"
          font-size="48" font-weight="820" fill="#222734" letter-spacing="0">NearRecurrences.jl</text>
    <text x="2" y="90" font-family="Inter, Avenir, Helvetica, Arial, sans-serif"
          font-size="21" font-weight="700" fill="#4f596b" letter-spacing="0">Near recurrences of a simulation, as it runs,</text>
    <text x="2" y="118" font-family="Inter, Avenir, Helvetica, Arial, sans-serif"
          font-size="21" font-weight="700" fill="#4f596b" letter-spacing="0">with its own states as the seeds of a search.</text>
  </g>
</svg>
"""

write(joinpath(@__DIR__, "logo.svg"), svg)
@printf("written assets/logo.svg, a stretch of period %.2f\n", 0.05*lag)
