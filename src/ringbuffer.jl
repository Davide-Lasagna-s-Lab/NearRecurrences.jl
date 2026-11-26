# The ring buffer, the last n samples of a stream, each indexed by its number in the stream.
#
# The buffer holds n copies of a template, made once. push! copies a new sample into the slot of the
# oldest, so that the stream allocates no state after the start. Sample k sits in slot mod1(k, n),
# and the samples held are firstindex(B) to lastindex(B), the last n pushed. With n = 4 after six
# samples
#
#     slot       1   2   3   4
#     sample     5   6   3   4      firstindex 3, lastindex 6
#
# and the next push overwrites sample 3 with sample 7. Reading a sample that has left the buffer, or
# that has not come yet, throws a BoundsError, so that an old sample is never read in place of
# another.

"""
    RingBuffer(x, n)

The last `n` samples of a stream, indexed by their number in the stream, from 1. It holds `n`
copies of `x`, made once. `push!(B, x)` copies `x` into the slot of the oldest sample, and `B[k]`
is sample `k`, the k-th pushed, while it is held, from `firstindex(B)` to `lastindex(B)`.

# Arguments
- `x`: a template of a sample, with `copy` and `copyto!`, as an array.
- `n::Int`: the number of samples held.
"""
mutable struct RingBuffer{X}
    slots::Vector{X} # the samples, sample k in slot mod1(k, n)
     last::Int       # the number of the last sample pushed, 0 before the first

    RingBuffer(x, n::Int) = new{typeof(x)}([copy(x) for _ in 1:n],
                                           0)
end

# ---- a new sample, copied into the slot of the oldest ----
function Base.push!(B::RingBuffer, x)
    B.last += 1
    copyto!(B.slots[mod1(B.last, length(B.slots))], x)
    return B
end

# ---- the samples held, the last n pushed ----
Base.firstindex(B::RingBuffer) = max(1, B.last - length(B.slots) + 1)
Base.lastindex(B::RingBuffer)  = B.last

# ---- sample k, while it is held ----
function Base.getindex(B::RingBuffer, k::Int)
    firstindex(B) <= k <= lastindex(B) ||
        throw(BoundsError(B, k))
    return B.slots[mod1(k, length(B.slots))]
end
