# Near recurrences of a simulation, searched as the simulation runs. The samples u₁, u₂, … come one
# every Δs. A near recurrence is a pair of samples, a start u(t) and an end u(t + T), whose distance
# d(t, T), minimised over the symmetries of the system by the user's function, is a local minimum in
# the plane (t, T), lies below a threshold, and is the smallest of the minima within a given apart,
# for periods T in given windows. The package keeps the last samples in a ring buffer, the last
# three rows d(i, ·) of each window, the candidates still undecided, and a histogram of the distances
# below the threshold, from which the probability of a distance follows (recurrences.jl). The states
# between the start and the end of a recurrence, its stretch, stay in the buffer, and the shooting
# points of a search for a periodic orbit are copies of them, the simulation's own states (seeds.jl).
module NearRecurrences

export Recurrences,
       seeds,
       finish!,
       probability

include("ringbuffer.jl")
include("recurrences.jl")
include("seeds.jl")

end
