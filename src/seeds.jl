# The seeds of a recurrence, the shooting points of a search for a periodic orbit by multiple
# shooting.
#
# A search with m segments needs m shooting points u₁, …, u_m along the stretch of the recurrence,
# from its start to its end. Let them be samples of the simulation, and let segment j advance its
# point by the same calls of the simulation that led from it to the next point. Each interior
# junction, the end of segment j minus u_{j+1}, then vanishes, and the closing junction, the end of
# segment m shifted minus u₁, is the end of the stretch shifted minus its start. The first residual
# of the search is then the recurrence itself and nothing else, of norm d ‖u₁‖ for the distance
# ‖S u_end - u_start‖/‖u_start‖. So the seeds are copies of samples held in the buffer, nothing
# interpolated or advanced, and every segment spans a whole number of samples.
#
# The offsets. Segment j goes from the sample start + o_j to the sample start + o_{j+1}, with
#
#     o_j = round((j - 1) ℓ/m),   j = 1, …, m + 1,
#
# and ℓ the lag, so that o₁ = 0 and o_{m+1} = ℓ, the quotient m ℓ/m being exact. Consecutive values
# (j - 1) ℓ/m differ by ℓ/m, and rounding moves each by half a sample at most. A segment thus spans
# a whole number of samples between ℓ/m - 1 and ℓ/m + 1, that is ℓ/m itself when it is whole, and
# its floor or its ceiling otherwise, one sample at least since m ≤ ℓ. With ℓ = 10 and m = 4 the
# offsets are 0, 2, 5, 8 and 10, the values 2.5 and 7.5 rounded to their even neighbours, and the
# segments have 2, 3, 3 and 2 samples. With m = ℓ every sample is a seed.

"""
    seeds(R, r, m)

The `m` shooting points of the recurrence `r` and the end of its stretch, copies of samples held in
the buffer, the simulation's own states bit for bit. Segment `j` goes from the sample `offsets[j]`
after the start to the sample `offsets[j + 1]`. Advancing `states[j]` by the calls of the simulation
over those samples gives `states[j + 1]`, and `distance(states[1], states[m + 1])` returns the
distance of `r` and its shifts, those that bring the end back onto the start. Nothing is
interpolated or advanced.

# Arguments
- `R::Recurrences`: the search that returned `r`.
- `r`: a recurrence of `R`, from `push!` or `finish!`, taken before the next `push!`. An
  `ArgumentError` is thrown once its stretch has left the buffer.
- `m::Int`: the number of segments, from 1 to `r.lag`. With `m == r.lag` every sample is a seed.

Returns `(states, offsets)`, with
- `states`: `m + 1` copies of samples, the shooting points `states[1:m]` and the end of the stretch
  `states[m + 1]`;
- `offsets`: their distances from the start, in samples, `round(Int, (j - 1)*r.lag/m)` for `j = 1`
  to `m + 1`, from `0` to `r.lag`.
"""
function seeds(R::Recurrences, r, m::Int)
    # ---- input checks ----
    1 <= m <= r.lag ||
        throw(ArgumentError("m must be from 1 to the lag of the recurrence, $(r.lag)"))
    firstindex(R.buffer) <= r.start && r.start + r.lag <= lastindex(R.buffer) ||
        throw(ArgumentError("the stretch of the recurrence has left the buffer, " *
                            "take the seeds before the next push!"))

    # ---- the offsets of the m + 1 states from the start, 0 first and the lag last ----
    offsets = [round(Int, (j - 1)*r.lag/m) for j in 1:(m + 1)]

    # ---- copies of the samples ----
    states = [copy(R.buffer[r.start + o]) for o in offsets]

    return (states  = states,
            offsets = offsets)
end
