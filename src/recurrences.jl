# The search of near recurrences as a stream, one sample at a time, with the whole trajectory never
# stored.
#
# The distance. distance(a, b) is that of a start a and an end b, minimised over the symmetries of
# the system, continuous or discrete, as many as the system has. It returns the distance alone, or
# a tuple (d, s₁, s₂, …) with the shifts that reach it, (d, σ, sy) for a shift σ along x and a
# translation sy along y. The search sees the distance only. A recurrence is a local minimum of the
# distance already minimised over the symmetries, and it carries the shifts of its pair, whatever
# the shifts of its neighbours.
#
# The rows. Sample k is the state at the time t_k = (k - 1) Δs. A window holds the periods from lo to
# hi, that is the lags ℓ from round(lo/Δs) to round(hi/Δs) in samples. With h the longest lag of all
# the windows, the sample k completes the row of the start i = k - h in every window,
#
#     d(i, ℓ) = distance(u_i, u_{i+ℓ}),   ℓ in the lags of the window,
#
# since every lag of every window is at most h. All the windows thus hold the rows of the same
# starts, and the candidates come start after start whatever their window.
#
# The candidates. Three rows are kept per window, those of i - 2, i - 1 and i. An entry d(i - 1, ℓ) of
# the middle row is a candidate if it is below the threshold of the window and below its eight
# neighbours, the entries of the three rows at ℓ - 1, ℓ and ℓ + 1 but itself. The first and the last
# lag of a window lack the neighbours on one side and are never candidates, nor are the first and the
# last row of the stream.
#
# The rule of the closest. With W = round(apart/Δs), two candidates are close if their starts differ
# by W samples at most and their lags too, in any windows. Two close minima share most of their
# stretch of trajectory, so they are repetitions of one guess. A candidate is kept unless a close one
# comes before it in the order of (distance, start, lag, window), so that of two close candidates the
# one of smaller distance is kept and the other is dropped, whatever their shifts. The kept ones are
# never close to one another. The last three keys break ties of distance. A tie arises where windows
# overlap and one pair is a candidate of both, which the rule keeps once, in the window of lower
# index. The rule depends on the candidates alone, not on the order in which they come.
#
# When a candidate is decided. A candidate of the start s depends on the candidates of the starts
# s - W to s + W. These have all come once the middle row has reached s + W + 1, since later
# candidates start at the middle row or after it. The sample k, whose middle row is k - h - 1, thus
# decides the candidates of the start s = k - h - 2 - W. push! returns those kept, and finish!
# decides the candidates still undecided when the stream stops. A dropped candidate stays among the
# candidates while a later decision needs it, that is while its start is above s - W.
#
# The length of the buffer. The recurrences decided at the sample k start at k - h - 2 - W, and their
# stretches, the samples from the start to the start plus the lag, at most k, must be held until the
# next push. The buffer thus holds the samples k - h - 2 - W to k, that is
#
#     N = h + W + 3
#
# samples. The next push writes sample k + 1 over sample k + 1 - N = k - h - 2 - W, the start of the
# oldest stretch. The test "Lorenz, two windows" checks that the recurrences start at the first
# sample held, and that the next push takes their stretch away. W may exceed the width of a window
# and its lags. It only delays the decisions and lengthens the buffer.
#
# The memory. The candidates held start within the 2W + 2 rows from s - W + 1 to the middle row.
# Being minima below the threshold, they are few, and their vector grows only when they outnumber
# all those held before. The recurrences kept at one start are not close, so their lags differ by more
# than W, and at most cld(h - ℓ_lo + 1, W + 1) of them are decided at one sample, with ℓ_lo the lowest
# lag. Their vector is given that capacity once. push! thus allocates nothing after the start, save
# what the invariant of the user returns and the growth of the candidates to a new largest number
# (test "a long stream").
#
# Windows that touch. Two windows (lo, T₁) and (T₁, hi) share the lag round(T₁/Δs), the last of one
# and the first of the other. A lag at the edge of a window lacks the neighbours on one side, so a
# minimum at the common lag is a candidate of neither window. Windows that overlap by two lags or
# more find it, and the rule of the closest keeps it once (test "a long stream").
#
# The invariant bound. An invariant v(u) is a real array that the symmetries leave unchanged,
# v(S u) = v(u), and that does not stretch distances, ‖v(a) - v(b)‖ ≤ ‖a - b‖, with ‖v(a)‖ = ‖a‖. For
# every symmetry S,
#
#     ‖v(b) - v(a)‖ = ‖v(S b) - v(a)‖ ≤ ‖S b - a‖,
#
# and so ‖v(b) - v(a)‖ is at most the minimum over S. Divided by ‖v(a)‖ = ‖a‖ this is the bound
# β(a, b) = ‖v(a) - v(b)‖/‖v(a)‖ ≤ d(a, b), which costs one pass over the entries of the invariants.
# A pair whose β is not below the threshold holds β in its entry, with no call of the distance. Such
# an entry is at or above the threshold, as its distance is. It is then no candidate, and it lies
# above any candidate among its neighbours, as its distance does, since a candidate is below the
# threshold. The candidates are therefore those of the distances alone.
#
# The probability of a distance. The fraction C(d) of the pairs closer than d, the correlation
# integral, sets the threshold of a search by a probability p, C(d_p) = p. The rows already hold the
# distances of every pair at the lags of the windows, so each window counts its pairs and keeps a
# histogram of log10 d for the distances below its threshold. Bin j holds
#
#     edges[j] ≤ log10 d < edges[j + 1],   edges[j] = -2 + (j - 1)/100,
#
# and the first bin the distances below 10^-1.99 as well, up to the bin that holds the threshold.
# With an invariant, a pair whose bound is at or above the threshold is not compared in full but
# counts as a pair. Every pair closer than the threshold has β ≤ d below it, so it is compared in
# full and binned, and the histogram is the same as without the invariant. For a d at an edge
# 10^-1.99, 10^-1.98, … at or below the smallest threshold, the pairs closer than d are those of the
# bins below it, in every window, and C(d) is exact. With the edges 0.01 apart in log10 d, a d
# between two edges lies 2.3 per cent at most above the one below, whose fraction it takes.

"""
    Recurrences(x, distance; spacing, windows, threshold, apart, invariant=nothing)

The search of near recurrences along a stream of samples, one every `spacing`. `push!(R, x)` takes
the next sample and returns the recurrences that it makes final, and `finish!(R)` those still
undecided when the stream stops. A recurrence is a pair of samples, a start and an end a period `T`
later in one of the windows, whose distance is a local minimum over the start and the period, lies
below the threshold of its window, and is the smallest of the minima within `apart` of it, across
the windows. The buffer holds the last samples, enough that the stretch of a recurrence, the
samples from its start to its end, is held until the next `push!`, for `seeds`.

# Arguments
- `x`: a template of a sample, with `copy` and `copyto!`, as an array. The buffer holds copies of
  it.
- `distance`: the distance of a start `a` and an end `b`, called as `distance(a, b)`, minimised over
  the symmetries of the system, continuous and discrete, as many as it has. It returns a number, or
  a tuple `(d, s₁, s₂, …)` whose first entry is the distance and the others the shifts that reach
  it, as `(d, σ, sy)`. It is called once here, on `(x, x)`, for the type of the shifts.
- `spacing::Real`: the time between two samples.
- `windows`: the periods searched, one tuple `(lo, hi)` per window, in time units. Its lags are
  `round(Int, lo/spacing)` to `round(Int, hi/spacing)` samples.
- `threshold`: the largest distance of a recurrence, excluded, one number for all the windows or one
  per window.
- `apart::Real`: two minima whose starts differ by `apart` at most and whose periods too, in time
  units, share most of their stretch and repeat one guess. Of two such minima, in any windows and
  with any shifts, the one of smaller distance is kept. Rounded to whole samples. It may exceed the
  width of a window.
- `invariant=nothing`: a function of a sample that returns a real array `v`, unchanged by the
  symmetries of `distance`, with `norm(v(a)) == norm(a)` and `norm(v(a) - v(b)) <= norm(a - b)` in the
  norm of `distance`. Then `norm(v(a) - v(b))/norm(v(a)) <= distance(a, b)`, and a pair whose bound
  is not below the threshold is not passed to `distance`. The recurrences are the same.

A recurrence is the named tuple `(t, T, d, shifts, window, start, lag)`, with the time of the start,
the period, the distance, the shifts that `distance` returns after it (an empty tuple if none), the
index of the window, and the start and the lag in samples. The k-th sample pushed is sample `k`, at
the time `t = (k - 1)*spacing`, and the period is `T = lag*spacing`. A recurrence is decided once
the longest window has moved `apart` past its start, so the recurrences of every window come in step
with the longest one.

Each window counts the pairs it compares, `R.pairs[w]`, and the distances below its threshold by
bins of `log10 d`, `R.counts[w][j]` for `R.edges[w][j] <= log10(d) < R.edges[w][j + 1]`, with the
edges `-2, -1.99, …` up to the bin of the threshold and the distances below `10^-1.99` in the first
bin. `probability(R, d)` reads them. The counts of several runs with the same thresholds add up bin
by bin, and so do their pairs.
"""
struct Recurrences{X, V, D, I, C}
        buffer::RingBuffer{X}                 # the last samples, sample k at index k
    invariants::V                             # their invariants in a RingBuffer, or nothing
      distance::D                             # distance(a, b), a number or a tuple (d, s₁, …)
     invariant::I                             # invariant(a), or nothing
       spacing::Float64                       # the time between two samples
          lags::Vector{UnitRange{Int}}        # the lags of each window, in samples
       longest::Int                           # the longest lag of all the windows, h
     threshold::Vector{Float64}               # the largest distance of a recurrence, per window
         apart::Int                           # the rule of the closest, W, in samples
          rows::Vector{Vector{Vector{Float64}}} # per window, the rows of the last three starts
    candidates::Vector{@NamedTuple{d::Float64, start::Int, lag::Int, window::Int}} # those held
         final::Vector{C}                     # the recurrences decided by the last sample
         pairs::Vector{Int}                   # per window, the pairs compared
         edges::Vector{Vector{Float64}}       # per window, the edges of the bins of log10 d
        counts::Vector{Vector{Int}}           # per window, the distances below the threshold by bin

    function Recurrences(x, distance;
                         spacing::Real,
                         windows,
                         threshold,
                         apart::Real,
                         invariant=nothing)
        # ---- input checks ----
        length(threshold) in (1, length(windows)) ||
            throw(ArgumentError("threshold must be one distance or one per window"))

        # ---- the windows in samples, their thresholds, and the rule of the closest ----
        lags   = [round(Int, lo/spacing):round(Int, hi/spacing) for (lo, hi) in windows]
        h      = maximum(last, lags)
        limits = length(threshold) == 1 ? fill(Float64(first(threshold)), length(lags)) :
                                          Float64.(collect(threshold))
        W      = round(Int, apart/spacing)

        # ---- the length of the buffer ----
        # The sample k decides the candidates of the start k - h - 2 - W, whose stretches reach back
        # to that sample, so the buffer holds h + W + 3 samples (see the top of this file).
        n = h + W + 3

        # ---- the samples, and their invariants if given ----
        buffer     = RingBuffer(x, n)
        invariants = invariant === nothing ? nothing : RingBuffer(invariant(x), n)

        # ---- three rows per window, those of the last three starts with the oldest first, at Inf
        #      until the stream fills them ----
        rows = [[fill(Inf, length(l)) for _ in 1:3] for l in lags]

        # ---- the type of a recurrence, with the shifts that distance returns after the distance ----
        v = distance(x, x)
        S = typeof(v isa Tuple ? Base.tail(v) : ())
        C = @NamedTuple{t::Float64,
                        T::Float64,
                        d::Float64,
                        shifts::S,
                        window::Int,
                        start::Int,
                        lag::Int}

        # ---- the candidates, and the recurrences decided by one sample, cld(h - ℓ_lo + 1, W + 1) at
        #      most (see the top of this file) ----
        candidates = @NamedTuple{d::Float64, start::Int, lag::Int, window::Int}[]
        final      = C[]
        sizehint!(final, cld(h - minimum(first, lags) + 1, W + 1))

        # ---- the histograms of log10 d, bins of 0.01 from -2 to the bin of the threshold, which
        #      bins = ceil(100 (log10 threshold + 2)) of them reach (see the top of this file) ----
        bins   = [max(1, ceil(Int, 100*(log10(limit) + 2))) for limit in limits]
        edges  = [[(j - 200)/100 for j in 0:J] for J in bins]
        counts = [zeros(Int, J) for J in bins]
        pairs  = zeros(Int, length(lags))

        return new{typeof(x),
                   typeof(invariants),
                   typeof(distance),
                   typeof(invariant),
                   C}(buffer,
                      invariants,
                      distance,
                      invariant,
                      Float64(spacing),
                      lags,
                      h,
                      limits,
                      W,
                      rows,
                      candidates,
                      final,
                      pairs,
                      edges,
                      counts)
    end
end

"""
    push!(R::Recurrences, x)

Takes `x` as the next sample of the stream and returns the recurrences that it makes final, a
vector that is empty most of the time. `x` is copied into the buffer. The sample completes the row
of one start in every window, the distances from it to the samples at the lags of the window, the
same start for all the windows, set by the longest lag. The middle of the last three rows of each
window is then tested for local minima below the threshold, which join the candidates. The sample
also decides the candidates that start `apart` and one more sample before the middle row, since
no later candidate can be close to them, and returns those kept by the rule of the closest. The stretch
of a recurrence returned is held until the next `push!`, for `seeds`. So is the vector returned, the
same at every call, which the next `push!` empties. Copy the recurrences to keep them, as
`append!(list, push!(R, x))` does. After the start `push!` allocates nothing, save the array that
`invariant` returns and the growth of the candidates to a new largest number.

# Arguments
- `R::Recurrences`: the search, changed.
- `x`: the next sample, not changed.
"""
function Base.push!(R::Recurrences, x)
    # ---- the sample and its invariant, into the buffers ----
    push!(R.buffer, x)
    R.invariant === nothing ||
        push!(R.invariants, R.invariant(x))
    k = lastindex(R.buffer)

    # ---- the recurrences that this sample makes final, in the vector of the last sample ----
    final = empty!(R.final)

    # ---- the start whose lags all have their sample in every window ----
    i = k - R.longest
    i >= 1 ||
        return final

    # ---- its row in every window, into the storage of the oldest row, so that the rows are those
    #      of i - 2, i - 1 and i ----
    for w in eachindex(R.lags)
        rows    = R.rows[w]
        oldest  = rows[1]
        rows[1] = rows[2]
        rows[2] = rows[3]
        rows[3] = oldest
        row!(rows[3], R, w, i)
        R.pairs[w] += length(R.lags[w])
    end

    # ---- the middle row i - 1 needs a row before it ----
    i >= 3 ||
        return final

    # ---- the candidates of the start s = i - 2 - W are decided, all those within W of them having
    #      come ----
    s = i - 2 - R.apart
    decide!(final, R, s:s)

    # ---- those that start at s - W or before leave, no later decision needing them, moved out in
    #      place, since filter! would give up the capacity of the vector ----
    kept = 0
    for c in R.candidates
        c.start > s - R.apart ||
            continue
        kept              += 1
        R.candidates[kept] = c
    end
    resize!(R.candidates, kept)

    # ---- the local minima of the middle row of every window join the candidates ----
    for w in eachindex(R.lags)
        minima!(R, w, i - 1)
    end

    return final
end

"""
    finish!(R::Recurrences)

The recurrences still undecided when the stream stops, a vector in the order of their starts. No
sample follows, so every candidate held is decided by the rule of the closest among those held.
Their stretches are held until the next `push!`, for `seeds`. The candidates are then emptied, so
that a second call returns none.

# Arguments
- `R::Recurrences`: the search, changed.
"""
function finish!(R::Recurrences)
    # ---- the starts decided so far, up to that of the last sample ----
    decided = lastindex(R.buffer) - R.longest - 2 - R.apart

    # ---- the others, all of them ----
    rest = empty(R.final)
    decide!(rest, R, (decided + 1):typemax(Int))
    empty!(R.candidates)

    return rest
end

"""
    probability(R::Recurrences, d)

The fraction of the pairs compared so far that are closer than `d`, pooled over the windows, the
correlation integral of the pairs at the lags searched. It reads the histograms of the windows.
When `log10(d)` is an edge of their bins, `-1.99, -1.98, …`, the fraction is exact, since every
pair closer than the threshold is compared in full, with the invariant or without. Between two
edges it is the fraction closer than the edge below `d`, at most 2.3 per cent smaller in distance.
A pair at a lag of two windows counts in both.

# Arguments
- `R::Recurrences`: the search.
- `d::Real`: the distance, at most the smallest threshold, below which the histograms are complete.
"""
function probability(R::Recurrences, d::Real)
    # ---- input checks ----
    d <= minimum(R.threshold) ||
        throw(ArgumentError("d must be at most the smallest threshold, $(minimum(R.threshold)), " *
                            "the histograms holding only the distances below it"))

    # ---- the pairs in the bins below the one of d, which lie below its edge edges[j] ≤ log10 d ----
    closer = 0
    for w in eachindex(R.counts)
        j       = searchsortedlast(R.edges[w], log10(d))
        closer += sum(R.counts[w][1:(j - 1)])
    end

    return closer/sum(R.pairs)
end


# ---------------------------------------------------------------------------- #
# the steps of a sample, the row, the minima and the decision                  #
# ---------------------------------------------------------------------------- #

# The row of the start i in the window w, into row: the distances from sample i to the samples at the
# lags of w. With an invariant, the bound β of each pair first, by a loop over the entries of the
# invariants that allocates nothing, and the distance only where β is below the threshold. A
# distance below the threshold goes into its bin of the histogram of w.
function row!(row::Vector{Float64}, R::Recurrences, w::Int, i::Int)
    lags      = R.lags[w]
    threshold = R.threshold[w]
    a         = R.buffer[i]

    # ---- the invariant of the start and its norm, for the bounds ----
    va = R.invariant === nothing ? nothing : R.invariants[i]
    na = R.invariant === nothing ? 0.0 : sqrt(sum(abs2, va))

    for l in eachindex(lags)
        j = i + lags[l]

        # ---- the bound β = ‖v(a) - v(b)‖/‖v(a)‖, which the pair keeps if it is not below the
        #      threshold ----
        if R.invariant !== nothing
            vb = R.invariants[j]
            e  = 0.0
            for q in eachindex(va, vb)
                e += (va[q] - vb[q])^2
            end
            β = sqrt(e)/na

            if β >= threshold
                row[l] = β
                continue
            end
        end

        # ---- the distance in full, and its bin if it is below the threshold ----
        d      = first(R.distance(a, R.buffer[j]))
        row[l] = d
        if d < threshold
            k = clamp(searchsortedlast(R.edges[w], log10(d)), 1, length(R.counts[w]))
            R.counts[w][k] += 1
        end
    end

    return row
end

# The local minima of the middle row of the window w, that of the start i, below the threshold of w,
# into the candidates.
function minima!(R::Recurrences, w::Int, i::Int)
    before, middle, after = R.rows[w]
    lags = R.lags[w]

    for l in 2:(length(lags) - 1)
        d = middle[l]

        # ---- below the threshold, and below its eight neighbours ----
        d < R.threshold[w] ||
            continue
        d < min(before[l - 1], before[l], before[l + 1],
                middle[l - 1],            middle[l + 1],
                 after[l - 1],  after[l],  after[l + 1]) ||
            continue

        push!(R.candidates, (d      = d,
                             start  = i,
                             lag    = lags[l],
                             window = w))
    end

    return R.candidates
end

# The candidates whose start lies in starts, decided by the rule of the closest, into out: a
# candidate is kept unless a candidate within W in start and in lag comes before it in the order of
# (distance, start, lag, window). The candidates held are those that the decision needs. The shifts
# of a recurrence kept come from one more call of the distance on its pair, since the rows hold
# distances only.
function decide!(out::Vector, R::Recurrences, starts::UnitRange{Int})
    for c in R.candidates
        c.start in starts ||
            continue

        # ---- a close candidate before it in the order drops it ----
        before(b) = abs(b.start - c.start) <= R.apart &&
                    abs(b.lag - c.lag) <= R.apart &&
                    (b.d, b.start, b.lag, b.window) < (c.d, c.start, c.lag, c.window)
        any(before, R.candidates) &&
            continue

        # ---- kept, with its shifts ----
        v = R.distance(R.buffer[c.start], R.buffer[c.start + c.lag])
        push!(out, (t      = (c.start - 1)*R.spacing,
                    T      = c.lag*R.spacing,
                    d      = c.d,
                    shifts = v isa Tuple ? Base.tail(v) : (),
                    window = c.window,
                    start  = c.start,
                    lag    = c.lag))
    end

    return out
end
