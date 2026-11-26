# The tests: the ring buffer, then the stream against a reference built from the whole trajectory, on
# the Lorenz system and on the Kuramoto–Sivashinsky equation with its two symmetries, the buffer, the
# seeds, finish!, the invariant bound, the probability, and a long stream for the runs of days.
#
#     julia --project=. -e 'using Pkg; Pkg.test()'
#
# The heavy cases, the Kuramoto–Sivashinsky stream of 10⁴ samples and the long stream of 5 × 10⁵,
# run only when the environment asks for them, NEARRECURRENCES_LONG=1. By default the same checks run
# on shorter streams, 4000 and 2 × 10⁴ samples.
import NearRecurrences
import LinearAlgebra
using Test

include(joinpath(@__DIR__, "systems.jl"))

const LONG = get(ENV, "NEARRECURRENCES_LONG", "0") == "1"


# ---------------------------------------------------------------------------- #
# the reference, from the whole trajectory                                     #
# ---------------------------------------------------------------------------- #

# The recurrences of a whole trajectory, with the rules written plainly. Every window has the rows of
# the starts whose lags all lie in the trajectory in the longest window, 1 to n. A candidate is an
# entry below the threshold of its window and below its eight neighbours, the first and the last
# start and lag excluded. A candidate is a recurrence unless a candidate of any window within W
# samples in start and in lag comes before it in the order of (distance, start, lag, window). The
# matrices are built by blocks of starts, each with one more start on either side for the
# neighbours, so that a long trajectory never needs its whole matrix. The histogram bins the entries
# below the threshold by the edges at or below their log10, the first bin taking those below -1.99.
# Returns the recurrences, the candidates, the histograms and the pairs of the windows.
function reference(samples, distance; spacing, windows, threshold, apart, block=10_000)
    W      = round(Int, apart/spacing)
    limits = length(threshold) == 1 ? fill(first(threshold), length(windows)) : threshold
    lags   = [round(Int, lo/spacing):round(Int, hi/spacing) for (lo, hi) in windows]
    n      = length(samples) - maximum(last, lags)

    candidates = []
    counts     = Vector{Int}[]
    pairs      = Int[]

    for (w, ℓs) in enumerate(lags)
        edges = [(j - 200)/100 for j in 0:max(1, ceil(Int, 100*(log10(limits[w]) + 2)))]
        bins  = zeros(Int, length(edges) - 1)

        for from in 1:block:n
            to     = min(from + block - 1, n)
            starts = max(1, from - 1):min(n, to + 1)
            D      = [first(distance(samples[i], samples[i + ℓ])) for i in starts, ℓ in ℓs]
            row(i) = i - first(starts) + 1

            # ---- the histogram of the starts of the block ----
            for i in from:to, c in eachindex(ℓs)
                x = D[row(i), c]
                x < limits[w] ||
                    continue
                bins[clamp(count(e -> e <= log10(x), edges), 1, length(bins))] += 1
            end

            # ---- the candidates of the starts of the block ----
            for i in max(2, from):min(n - 1, to), c in 2:(length(ℓs) - 1)
                d          = D[row(i), c]
                neighbours = [D[row(i) + a, c + b] for a in -1:1 for b in -1:1 if (a, b) != (0, 0)]
                d < limits[w] && d < minimum(neighbours) ||
                    continue
                push!(candidates, (d = d, start = i, lag = ℓs[c], window = w))
            end
        end
        push!(counts, bins)
        push!(pairs, n*length(ℓs))
    end

    # ---- the rule of the closest, over the candidates sorted by start ----
    sort!(candidates; by=c -> c.start)
    starts = [c.start for c in candidates]
    order(c) = (c.d, c.start, c.lag, c.window)
    found  = []
    for c in candidates
        near = candidates[searchsortedfirst(starts, c.start - W):searchsortedlast(starts, c.start + W)]
        any(b -> abs(b.lag - c.lag) <= W && order(b) < order(c), near) &&
            continue

        v = distance(samples[c.start], samples[c.start + c.lag])
        push!(found, (t      = (c.start - 1)*spacing,
                      T      = c.lag*spacing,
                      d      = c.d,
                      shifts = v isa Tuple ? Base.tail(v) : (),
                      window = c.window,
                      start  = c.start,
                      lag    = c.lag))
    end

    return found, candidates, counts, pairs
end

# ---- the recurrences in one order, by start, lag and window ----
ordered(list) = sort(list; by=r -> (r.start, r.lag, r.window))

# The seeds of r with m segments, while its stretch is held. Each state is its sample bit for bit,
# the calls of advance! over a segment lead from one state to the next bit for bit, and the first and
# the last state give the distance of r and its shifts exactly.
function check_seeds(R, r, m, samples, advance!, distance)
    g = NearRecurrences.seeds(R, r, m)

    # ---- the offsets ----
    @test length(g.states) == m + 1
    @test g.offsets == [round(Int, (j - 1)*r.lag/m) for j in 1:(m + 1)]
    @test first(g.offsets) == 0
    @test last(g.offsets) == r.lag

    # ---- the samples themselves, as copies ----
    @test all(g.states[j] == samples[r.start + g.offsets[j]] for j in 1:(m + 1))
    @test g.states[1] !== R.buffer[r.start]

    # ---- each segment repeats the calls of the simulation ----
    advanced = falses(m)
    for j in 1:m
        u = copy(g.states[j])
        for _ in 1:(g.offsets[j + 1] - g.offsets[j])
            advance!(u)
        end
        advanced[j] = u == g.states[j + 1]
    end
    @test all(advanced)

    # ---- the closing is the recurrence, with its shifts ----
    v = distance(g.states[1], g.states[m + 1])
    @test first(v) == r.d
    @test (v isa Tuple ? Base.tail(v) : ()) == r.shifts

    # ---- m out of range ----
    @test_throws "m must be" NearRecurrences.seeds(R, r, 0)
    @test_throws "m must be" NearRecurrences.seeds(R, r, r.lag + 1)
end

# The histograms and the probability of R against those of the reference. At the edges from -1.99
# to the smallest threshold, and at d between two edges, the probability is the fraction of the
# pairs in the bins below the edge at or under log10 d, in every window.
function check_probability(R, counts, pairs)
    @test R.counts == counts
    @test R.pairs == pairs

    smallest = minimum(R.threshold)
    for d in [exp10.(-1.99:0.01:log10(smallest)); smallest; range(0.0103, smallest; length=7)]
        closer = sum(eachindex(counts)) do w
            j = count(e -> e <= log10(d), R.edges[w])
            sum(counts[w][1:(j - 1)])
        end
        @test NearRecurrences.probability(R, d) == closer/sum(pairs)
    end
    @test_throws ArgumentError NearRecurrences.probability(R, 1.01*smallest)
end

# The live stream of n samples from x through R. Returns the recurrences in the order they come, the
# seeds of each with five segments, the samples, and for every sample after the tenth the bytes that
# push! allocates with the number of candidates it leaves.
function live!(R, advance!, x, n)
    samples    = Vector{Vector{Float64}}(undef, n)
    found      = []
    seeded     = []
    bytes      = zeros(Int, n)
    candidates = zeros(Int, n)

    for k in 1:n
        advance!(x)
        samples[k] = copy(x)

        # ---- the sample, the bytes that push! allocates, the candidates it leaves ----
        bytes[k]      = @allocated returned = push!(R, x)
        candidates[k] = length(R.candidates)

        # ---- the recurrences and their seeds, while their stretches are held ----
        for r in returned
            push!(found, r)
            push!(seeded, NearRecurrences.seeds(R, r, min(5, r.lag)))
        end
    end

    return found, seeded, samples, bytes, candidates
end

# ---- the same stream, from the samples stored: the recurrences and their seeds ----
function replay!(R, samples)
    found  = []
    seeded = []
    for x in samples
        for r in push!(R, x)
            push!(found, r)
            push!(seeded, NearRecurrences.seeds(R, r, min(5, r.lag)))
        end
    end
    return found, seeded
end

# ---- the seconds of a stream of the samples stored through R ----
function seconds!(R, samples)
    t = time_ns()
    for x in samples
        push!(R, x)
    end
    return (time_ns() - t)/1e9
end

# ---- the bytes allocated by the row of the start i in the window w, after a call that compiles ----
function row_allocations(R, w, i)
    row = zeros(length(R.lags[w]))
    NearRecurrences.row!(row, R, w, i)
    return @allocated NearRecurrences.row!(row, R, w, i)
end

# ---- two candidates within W samples in start and in lag ----
near(a, b, W) = abs(a.start - b.start) <= W && abs(a.lag - b.lag) <= W


# ---------------------------------------------------------------------------- #
# the tests                                                                    #
# ---------------------------------------------------------------------------- #

@testset "the ring buffer                          " begin
    B = NearRecurrences.RingBuffer([0.0], 4)
    @test firstindex(B) == 1
    @test lastindex(B) == 0

    # ---- ten samples through four slots, the last four held ----
    x = [0.0]
    for k in 1:10
        x[1] = k
        push!(B, x)
        @test firstindex(B) == max(1, k - 3)
        @test lastindex(B) == k
        @test all(B[j] == [j] for j in firstindex(B):lastindex(B))
        @test_throws BoundsError B[k + 1]
        @test_throws BoundsError B[firstindex(B) - 1]
    end

    # ---- after the wrap, sample k in slot mod1(k, 4) ----
    @test B.slots == [[9.0], [10.0], [7.0], [8.0]]

    # ---- push! copies, so that a change of x leaves the buffer ----
    x[1] = -1
    @test B[10] == [10.0]
end

@testset "input checks and the windows in samples    " begin
    # ---- the lags of a window are its periods over the spacing, rounded, both ends included, and
    #      apart in samples the same way: 0.3/0.1 and 0.7/0.1 fall just below 3 and 7 ----
    R = NearRecurrences.Recurrences([1.0, 1.0, 1.0], lorenz_distance;
                                    spacing   = 0.1,
                                    windows   = [(0.3, 0.7), (1.1, 1.7)],
                                    threshold = 0.1,
                                    apart     = 0.3)
    @test 0.3/0.1 < 3
    @test 0.7/0.1 < 7
    @test R.lags == [3:7, 11:17]
    @test R.apart == 3
    @test R.longest == 17
    @test length(R.buffer.slots) == 17 + 3 + 3

    @test_throws ArgumentError NearRecurrences.Recurrences([1.0, 1.0, 1.0], lorenz_distance;
                                                           spacing   = 0.05,
                                                           windows   = [(1.3, 1.8), (2.0, 2.6)],
                                                           threshold = [0.1, 0.1, 0.1],
                                                           apart     = 1.0)
end

@testset "Lorenz, two windows                      " begin
    # ---- apart of 20 samples, wider than the first window of 11 lags, a threshold per window ----
    settings = (spacing   = 0.05,
                windows   = [(1.3, 1.8), (2.0, 2.6)],
                threshold = [0.1, 0.15],
                apart     = 1.0)
    x = transient(lorenz_advance!, [1.0, 1.0, 1.0], 400)
    R = NearRecurrences.Recurrences(x, lorenz_distance; settings...)

    # ---- the stream, with the samples kept for the reference ----
    samples = Vector{Float64}[]
    found   = []
    oldest  = nothing # a recurrence whose start leaves the buffer at the next push
    for k in 1:20000
        lorenz_advance!(x)
        push!(samples, copy(x))
        final = push!(R, x)

        # ---- the buffer: the stretch that started at its first sample has left ----
        if oldest !== nothing
            @test_throws "left the buffer" NearRecurrences.seeds(R, oldest, 1)
            oldest = nothing
        end

        for r in final
            push!(found, r)

            # ---- the buffer: every recurrence starts at its first sample, in any window ----
            @test r.start == firstindex(R.buffer)
            oldest = r

            # ---- the seeds ----
            for m in (1, 3, 7, r.lag)
                check_seeds(R, r, m, samples, lorenz_advance!, lorenz_distance)
            end
        end
    end

    # ---- the recurrences still undecided, held for their seeds ----
    rest = NearRecurrences.finish!(R)
    for r in rest
        check_seeds(R, r, r.lag, samples, lorenz_advance!, lorenz_distance)
    end
    append!(found, rest)
    @test isempty(NearRecurrences.finish!(R))

    # ---- the same recurrences as the reference, in both windows, with the rule of the closest at
    #      work within a window and across the two ----
    expected, candidates, counts, pairs = reference(samples, lorenz_distance; settings...)
    W = R.apart
    @test length(found) == length(expected)
    @test ordered(found) == ordered(expected)
    @test all(w -> count(r -> r.window == w, found) >= 10, 1:2)
    @test length(candidates) > length(expected)
    @test any(c -> c.window == 2 && any(b -> b.window == 1 && near(b, c, W), candidates), candidates)
    @test all(r -> r.shifts == (), found)

    # ---- the histograms and the probability ----
    check_probability(R, counts, pairs)
end

@testset "Kuramoto–Sivashinsky, a shift and a reflection" begin
    # ---- apart of 40 samples, wider than both windows, a threshold per window ----
    settings = (spacing   = 0.25,
                windows   = [(9.0, 12.0), (30.0, 35.0)],
                threshold = [0.2, 0.15],
                apart     = 10.0)
    u = zeros(30)
    u[1] = 0.3
    u[4] = -0.2
    u[7] = 0.1
    u = transient(ks_advance!, u, 400)

    # ---- the test system: the distance is that of the shift it returns, and the invariant bounds it ----
    a          = copy(u)
    b          = transient(ks_advance!, copy(u), 47)
    d, σ, r    = ks_distance(a, b)
    @test LinearAlgebra.norm(ks_shift(b, σ, r) - a)/LinearAlgebra.norm(a) ≈ d
    bounded = map(1:200) do _
        ks_advance!(b)
        va = ks_invariant(a)
        vb = ks_invariant(b)
        LinearAlgebra.norm(va - vb) <= first(ks_distance(a, b))*LinearAlgebra.norm(va) + 1e-12
    end
    @test all(bounded)

    # ---- the stream with every distance in full, and with the invariant, on the same samples, with
    #      the calls of the distance counted ----
    calls = zeros(Int, 2)
    full  = NearRecurrences.Recurrences(u, (a, b) -> (calls[1] += 1; ks_distance(a, b));
                                        settings...)
    bound = NearRecurrences.Recurrences(u, (a, b) -> (calls[2] += 1; ks_distance(a, b));
                                        settings...,
                                        invariant = ks_invariant)
    samples = Vector{Float64}[]
    found   = ([], [])
    for k in 1:(LONG ? 10_000 : 4000)
        ks_advance!(u)
        push!(samples, copy(u))
        for (R, list) in zip((full, bound), found)
            for r in push!(R, u)
                push!(list, r)
                for m in (1, 4, r.lag)
                    check_seeds(R, r, m, samples, ks_advance!, ks_distance)
                end
            end
        end
    end
    append!(found[1], NearRecurrences.finish!(full))
    append!(found[2], NearRecurrences.finish!(bound))

    # ---- the same recurrences as the reference, with both shifts, with the invariant or without ----
    expected, candidates, counts, pairs = reference(samples, ks_distance; settings...)
    @test ordered(found[1]) == ordered(expected)
    @test ordered(found[2]) == ordered(expected)
    @test found[2] == found[1]
    @test all(w -> count(r -> r.window == w, found[1]) >= 3, 1:2)
    @test length(candidates) > length(expected)
    @test all(r -> length(r.shifts) == 2, found[1])
    @test Set(r.shifts[2] for r in found[1]) == Set([0, 1])

    # ---- the shifts close the stretch: S(σ) Rʳ of the end lies at d from the start ----
    @test all(found[1]) do r
        a = samples[r.start]
        b = samples[r.start + r.lag]
        LinearAlgebra.norm(ks_shift(b, r.shifts...) - a)/LinearAlgebra.norm(a) ≈ r.d
    end

    # ---- the histograms and the probability, with the invariant or without ----
    check_probability(full, counts, pairs)
    check_probability(bound, counts, pairs)

    # ---- the invariant spares most distances, and its row allocates nothing ----
    @test calls[2] < calls[1]/5
    @test row_allocations(bound, 1, lastindex(bound.buffer) - bound.longest) == 0

    # ---- the cost per pair, with every distance in full and with the invariant, on the samples
    #      stored, after a first stream that compiles ----
    cost = map((nothing, ks_invariant)) do invariant
        stream() = NearRecurrences.Recurrences(samples[1], ks_distance; settings...,
                                               invariant = invariant)
        seconds!(stream(), samples[1:500])
        R = stream()
        s = seconds!(R, samples[1:2000])
        s/sum(R.pairs)
    end
    println("Kuramoto–Sivashinsky: $(round(1e6*cost[1]; digits=2)) μs per pair with every distance ",
            "in full, $(round(1e6*cost[2]; digits=2)) μs with the invariant, which compares ",
            "$(round(calls[2]/calls[1]; sigdigits=2)) of the pairs in full")
end

@testset "a long stream                            " begin
    # ---- five windows of 11 to 201 lags, a threshold each, the first two touching at T = 1.8, the
    #      lag 36, the second and the third overlapping on the lags 44 to 52, the rotation of the
    #      Lorenz system, and apart of 100 samples, wider than all the windows but the last ----
    settings = (spacing   = 0.05,
                windows   = [(1.3, 1.8), (1.8, 2.6), (2.2, 3.0), (10.0, 12.0), (40.0, 50.0)],
                threshold = [0.1, 0.12, 0.12, 0.15, 0.2],
                apart     = 5.0)
    x = transient(lorenz_advance!, [1.0, 1.0, 1.0], 400)
    R = NearRecurrences.Recurrences(copy(x), rotated_distance; settings...)
    n = LONG ? 5*10^5 : 2*10^4

    # ---- the stream, live, across many wraps of the buffer ----
    found, seeded, samples, bytes, held = live!(R, lorenz_advance!, x, n)
    rest = NearRecurrences.finish!(R)
    @test n/length(R.buffer.slots) > 15

    # ---- the memory: push! allocates only when the candidates outnumber all before, and the
    #      recurrences decided by one sample stay within their capacity ----
    most = accumulate(max, held)
    grew = [held[k] > most[k - 1] for k in 11:n]
    @test all((bytes[11:n] .== 0) .| grew)
    @test sum(bytes[11:n] .> 0) <= 12
    @test length(R.final) <= cld(R.longest - minimum(first, R.lags) + 1, R.apart + 1)

    # ---- the same recurrences, histograms and pairs as the reference by blocks ----
    expected, candidates, counts, pairs = reference(samples, rotated_distance; settings...)
    W     = R.apart
    every = [found; rest]
    @test ordered(every) == ordered(expected)
    @test R.counts == counts
    @test R.pairs == pairs
    @test count(r -> r.window in (1, 2, 4, 5), every) >= 20

    # ---- the rotation: both shifts among the recurrences, and the rule of the closest drops
    #      candidates in favour of close ones of the other shift ----
    @test Set(r.shifts[1] for r in every) == Set([0, 1])
    rotations(c) = rotated_distance(samples[c.start], samples[c.start + c.lag])[2]
    @test any(candidates) do c
        any(b -> near(b, c, W) && b.d < c.d && rotations(b) != rotations(c), candidates)
    end

    # ---- the touching windows: the common lag is the edge of both, a candidate of neither, while
    #      one window over both has candidates there ----
    common = round(Int, 1.8/0.05)
    @test last(R.lags[1]) == common == first(R.lags[2])
    @test !any(c -> c.window in (1, 2) && c.lag == common, candidates)
    merged = reference(samples, rotated_distance; settings..., windows = [(1.3, 2.6)],
                      threshold = 0.12)[2]
    @test any(c -> c.lag == common, merged)

    # ---- the overlapping windows: a pair that is a candidate of both is one recurrence ----
    @test any(c -> c.window == 3 && any(b -> b.window == 2 && b.start == c.start && b.lag == c.lag,
                                        candidates), candidates)
    @test length(unique((r.start, r.lag) for r in every)) == length(every)

    # ---- the replay of the samples stored gives the same stream: the recurrences in the same
    #      order, their seeds bit for bit, and the state of the search ----
    again           = NearRecurrences.Recurrences(copy(samples[1]), rotated_distance; settings...)
    found₂, seeded₂ = replay!(again, samples)
    @test found₂ == found
    @test all(g.states == h.states && g.offsets == h.offsets for (g, h) in zip(seeded, seeded₂))
    @test again.buffer.slots == R.buffer.slots
    @test again.rows == R.rows
    @test again.counts == R.counts
    @test NearRecurrences.finish!(again) == rest
end

@testset "the seeds and the steps of the integrator" begin
    # A sample is advanced by one call over the span `spacing`, with steps of dt = 0.005 and a last
    # step that ends the span. A seed is advanced over a segment of n samples in three ways: by n calls
    # as the simulation, by one call over n*spacing with spacing = 10 dt, a whole number of steps, and
    # by one call over n*spacing with spacing = 10.5 dt, where the short steps fall elsewhere.
    dt = 0.005
    for (spacing, whole) in ((10dt, true), (10.5dt, false)), distance in (lorenz_distance, rotated_distance)
        advance! = x -> lorenz_span!(x, spacing, dt)
        x = transient(advance!, [1.0, 1.0, 1.0], 400)
        R = NearRecurrences.Recurrences(x, distance;
                                        spacing   = spacing,
                                        windows   = [(1.3, 1.8)],
                                        threshold = 0.15,
                                        apart     = 1.0)
        calls = Bool[]
        span  = Float64[]
        closing = Float64[]
        for k in 1:6000
            advance!(x)
            for r in push!(R, x)
                g = NearRecurrences.seeds(R, r, 4)
                for j in 1:4
                    n = g.offsets[j + 1] - g.offsets[j]

                    # ---- n calls of the simulation, the next seed bit for bit ----
                    u = copy(g.states[j])
                    for _ in 1:n
                        advance!(u)
                    end
                    push!(calls, u == g.states[j + 1])

                    # ---- one call over n*spacing, the next seed to rounding or not ----
                    v = lorenz_span!(copy(g.states[j]), n*spacing, dt)
                    push!(span, LinearAlgebra.norm(v - g.states[j + 1])/LinearAlgebra.norm(g.states[j + 1]))

                    # ---- the closing junction, from the last seed advanced, with the shift of r ----
                    if j == 4
                        e = distance(g.states[1], u)
                        @test first(e) == r.d
                        @test (e isa Tuple ? Base.tail(e) : ()) == r.shifts
                        b = distance === rotated_distance ? rotate(v, r.shifts[1]) : v
                        push!(closing, abs(lorenz_distance(g.states[1], b) - r.d))
                    end
                end
            end
        end
        @test length(calls) >= 20
        @test all(calls)
        println("spacing of $(spacing/dt) steps, $(nameof(distance)): $(length(calls)) segments, ",
                "one call over the segment misses the next seed by $(round(minimum(span); sigdigits=2)) ",
                "to $(round(maximum(span); sigdigits=2)), the closing by up to ",
                "$(round(maximum(closing); sigdigits=2))")
        if whole
            @test maximum(span) < 1e-12
            @test maximum(closing) < 1e-12
        else
            @test minimum(span) > 1e-12
        end
    end
end
