# The invariant bound on Kolmogorov flow, in the set-up of the searches of long orbits: Re = 40, the
# forcing wavenumber 4, the modes up to 33 on a grid of 100², the scheme CB3R2R3e with steps of 0.01,
# one field on the CPU, a sample every 0.1, from a smooth random field after a transient of 300 time
# units. The distance is that of BatchedKolmogorovFlow.jl, minimised over the shift along x and the
# translations along y, and the invariant its weighted moduli. A stream of 2000 time units through
# the search, for the windows 19 to 21 and 95 to 105 below 0.118 with apart 30, counts the pairs it
# compares in full. Alongside, for the figure only, the distance and the bound of the pairs of every
# fifth start of the first 500 time units, at every lag of the windows, without the skip, and the time
# of each on one core. Writes them to bound.csv, which figure.jl draws in an environment with
# CairoMakie. The two run apart, since CairoMakie brings a second FFT package whose plans make those of
# FFTW ambiguous in the propagator.
#
#     JULIA_LOAD_PATH="@:path/to/NearRecurrences.jl:@stdlib" \
#         julia --project=path/to/BatchedKolmogorovFlow.jl/examples examples/kolmogorov/bound.jl
#     julia --project=examples examples/kolmogorov/figure.jl
import NearRecurrences
import BatchedKolmogorovFlow
import Random
using Printf


# ---------------------------------------------------------------------------- #
# the stream and the pairs                                                     #
# ---------------------------------------------------------------------------- #

const n  = 33     # the largest wavenumber
const m  = 49     # the grid has 2m + 2 = 100 points along each direction
const Re = 40.0   # the Reynolds number
const kf = 4      # the wavenumber of the forcing

# ---- the bound of the pair (a, b) from their moduli, one pass over their entries ----
function bound(va, vb)
    e  = 0.0
    na = 0.0
    for q in eachindex(va, vb)
        e  += (va[q] - vb[q])^2
        na += va[q]^2
    end
    return sqrt(e/na)
end

# The stream of 2000 time units, the calls of the distance counted. Returns the recurrences, the
# search, the calls, the samples of the first 500 time units for the figure, and the distance.
function stream()
    G = BatchedKolmogorovFlow.NonlinearPropagator(n, m, Re, kf; B=1, dt=0.01)
    D = BatchedKolmogorovFlow.Distance(n, kf)

    # ---- a smooth random field, after a transient ----
    Random.seed!(1)
    x  = Random.randn!(BatchedKolmogorovFlow.FTFields(n, 1))
    x .= x ./ (1 .+ x.j.^2 .+ x.k.^2)
    BatchedKolmogorovFlow.apply_symmetry!(x)
    G(x, 300.0)

    # ---- the search, with the calls of the distance counted ----
    calls = Ref(0)
    R     = NearRecurrences.Recurrences(x, (a, b) -> (calls[] += 1; D(a, b));
                                        spacing   = 0.1,
                                        windows   = [(19.0, 21.0), (95.0, 105.0)],
                                        threshold = 0.118,
                                        apart     = 30.0,
                                        invariant = BatchedKolmogorovFlow.moduli)

    found   = []
    samples = typeof(x)[]
    for k in 1:20_000
        G(x, 0.1)                                     # one sample, ten steps
        k <= 5000 &&
            push!(samples, copy(x))
        append!(found, push!(R, x))
    end
    append!(found, NearRecurrences.finish!(R))

    return found, R, calls[], samples, D
end

found, R, calls, samples, D = stream()
pairs = sum(R.pairs)
@printf("%d recurrences in 2000 time units, %d pairs, %d compared in full, %.2f per cent\n",
        length(found), pairs, calls, 100calls/pairs)
for r in found
    @printf("  t = %7.1f, T = %6.1f, d = %.3f, σ = %+.3f, sy = %+d\n", r.t, r.T, r.d, r.shifts...)
end

# The distance and the bound of the pairs of every fifth start of the samples, at every lag of the
# windows, and the seconds that each takes in all
function sampled(samples, D)
    moduli = BatchedKolmogorovFlow.moduli.(samples)
    lags   = [190:210; 950:1050]
    d      = Float64[]
    β      = Float64[]
    td     = 0.0
    tβ     = 0.0
    for i in 1:5:(length(samples) - last(lags)), ℓ in lags
        t₀ = time_ns()
        push!(d, first(D(samples[i], samples[i + ℓ])))
        t₁ = time_ns()
        push!(β, bound(moduli[i], moduli[i + ℓ]))
        t₂ = time_ns()
        td += (t₁ - t₀)/1e9
        tβ += (t₂ - t₁)/1e9
    end
    return d, β, td, tβ
end

d, β, td, tβ = sampled(samples, D)
@printf("%d pairs: a distance %.1f μs, a bound %.2f μs, largest β/d %.4f\n",
        length(d), 1e6*td/length(d), 1e6*tβ/length(d), maximum(β ./ d))


# ---- the pairs, with the counts of the stream in the first line ----
open(joinpath(@__DIR__, "bound.csv"), "w") do io
    println(io, "# pairs $(pairs), compared $(calls), recurrences $(length(found)), ",
            "distance $(round(1e6*td/length(d); digits=2)) μs, bound $(round(1e6*tβ/length(d); digits=3)) μs")
    println(io, "distance,bound")
    for (x, y) in zip(d, β)
        println(io, x, ",", y)
    end
end
println("written examples/kolmogorov/bound.csv")
