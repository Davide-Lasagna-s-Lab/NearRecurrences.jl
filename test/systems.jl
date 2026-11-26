# The systems of the tests, each sampled by fixed steps of a Runge–Kutta scheme written here, so that a
# sample advanced by the same calls gives the next one bit for bit.
#
# The Lorenz system at σ = 10, ρ = 28, β = 8/3, sampled every 0.05 by five steps of 0.01 of the
# classical scheme. Two distances of a start a and an end b. The plain one is ‖b - a‖/‖a‖, a number.
# The other is minimised over the symmetry of the system, the rotation R(x, y, z) = (-x, -y, z) by π
# about the z axis, and returns (d, r) with r = 0 or 1 the rotations of b that reach d. Its
# invariant is v = (√(x² + y²), z), unchanged by R, with ‖v(a)‖ = ‖a‖. By the reverse triangle
# inequality |√(x_b² + y_b²) - √(x_a² + y_a²)| ≤ ‖(x_b, y_b) - (x_a, y_a)‖, so that
# ‖v(b) - v(a)‖ ≤ ‖b - a‖, and the same for R b in place of b, since v(R b) = v(b). The system is also
# advanced over a time span by steps of a nominal dt and a last step that ends the span, as the
# stepper of Flows.jl does, for the test of the seeds against the steps of the integrator.
#
# The Kuramoto–Sivashinsky equation u_t = -u u_x - u_xx - u_xxxx on a periodic domain of length
# L = 22, with the modes a_k, k = 1, …, 15, of u = Σ a_k e^{i q_k x} + c.c., q_k = 2πk/L, and the mean
# zero. The state is the real vector (Re a₁, Im a₁, …, Re a₁₅, Im a₁₅). The product u² is formed on
# 48 points, enough that no mode up to 30 folds onto the 15 kept, and the scheme is the classical
# one on the equation of e^{-Λt} a, the linear part Λ = q² - q⁴ taken exactly, with steps of 0.025,
# sampled every 0.25. The equation has two symmetries, the shift S(σ) u(x) = u(x + σ), which takes
# a_k to e^{i q_k σ} a_k, continuous, and the reflection R u(x) = -u(-x), which takes a_k to
# -conj(a_k), discrete. The distance is minimised over both and returns (d, σ, r). With the norm
# ‖a‖² = 2 Σ |a_k|²,
#
#     ‖S(σ) Rʳ b - a‖² = ‖a‖² + ‖b‖² - 2 f_r(σ),   f_r(σ) = 2 Re Σ_k conj(a_k) (Rʳ b)_k e^{i q_k σ},
#
# a trigonometric polynomial in σ for each r. It is taken at 128 shifts, and four Newton steps refine
# the largest over the shifts and r, as for Kolmogorov flow. Its invariant is v_k = √2 |a_k|,
# unchanged by both symmetries, of norm ‖a‖, and ||b_k| - |a_k|| ≤ |b_k - a_k| mode by mode.


# ---- one step of the classical Runge–Kutta scheme of f!(dx, x), in place, with the workspace k ----
function rk4!(f!, x, dt, k)
    k₁, k₂, k₃, k₄, y = k
    f!(k₁, x)
    y .= x .+ (dt/2) .* k₁
    f!(k₂, y)
    y .= x .+ (dt/2) .* k₂
    f!(k₃, y)
    y .= x .+ dt .* k₃
    f!(k₄, y)
    x .+= (dt/6) .* (k₁ .+ 2 .* k₂ .+ 2 .* k₃ .+ k₄)
    return x
end

# ---- the state after a transient of n samples from u ----
function transient(advance!, u, n)
    for _ in 1:n
        advance!(u)
    end
    return u
end


# ---------------------------------------------------------------------------- #
# Lorenz                                                                       #
# ---------------------------------------------------------------------------- #

function lorenz!(dx, x)
    dx[1] = 10*(x[2] - x[1])
    dx[2] = x[1]*(28 - x[3]) - x[2]
    dx[3] = x[1]*x[2] - 8/3*x[3]
    return dx
end

const LORENZ_WORK = ntuple(_ -> zeros(3), 5)

# ---- one sample ahead, five steps of 0.01 ----
function lorenz_advance!(x)
    for _ in 1:5
        rk4!(lorenz!, x, 0.01, LORENZ_WORK)
    end
    return x
end

# ---- ‖b - a‖/‖a‖ ----
lorenz_distance(a, b) = sqrt((b[1] - a[1])^2 + (b[2] - a[2])^2 + (b[3] - a[3])^2)/
                        sqrt(a[1]^2 + a[2]^2 + a[3]^2)

# ---- the least of ‖b - a‖ and ‖R b - a‖ over ‖a‖, and the rotations r of b that reach it ----
function rotated_distance(a, b)
    d₀ = sqrt((b[1] - a[1])^2 + (b[2] - a[2])^2 + (b[3] - a[3])^2)
    d₁ = sqrt((-b[1] - a[1])^2 + (-b[2] - a[2])^2 + (b[3] - a[3])^2)
    na = sqrt(a[1]^2 + a[2]^2 + a[3]^2)
    return d₁ < d₀ ? (d₁/na, 1) : (d₀/na, 0)
end

# ---- (√(x² + y²), z), unchanged by R, of norm ‖a‖ ----
rotated_invariant(a) = [hypot(a[1], a[2]), a[3]]

# ---- R b ----
rotate(b, r) = r == 0 ? copy(b) : [-b[1], -b[2], b[3]]

# The Lorenz system advanced over the span τ by the classical scheme with the steps of Flows.jl: the
# steps start at the times of the range 0:dt:τ and have the nominal dt, but the last one is τ - t,
# so that the span ends on τ. When τ is not a whole number of steps the last step is shorter. When
# it is, the range may stop a few ulps short of τ, and a last step shorter than a millionth of dt is
# such rounding, merged into the step before, which then differs from dt by rounding alone.
function lorenz_span!(x, τ, dt)
    times = 0:dt:τ
    lossy = length(times) == 1 || abs(τ - last(times)) > 1e-6*dt
    steps = lossy ? length(times) : length(times) - 1
    for i in 1:steps
        h = i == steps ? τ - times[i] : step(times)
        rk4!(lorenz!, x, h, LORENZ_WORK)
    end
    return x
end


# ---------------------------------------------------------------------------- #
# Kuramoto–Sivashinsky                                                         #
# ---------------------------------------------------------------------------- #

const KS_L    = 22.0                                                   # the length of the domain
const KS_K    = 15                                                     # the modes kept
const KS_N    = 48                                                     # the points of the products
const KS_DT   = 0.025                                                  # the time step
const KS_Q    = [2π*k/KS_L for k in 1:KS_K]                            # the wavenumbers
const KS_E    = [cis(2π*k*n/KS_N) for k in 1:KS_K, n in 0:(KS_N - 1)]  # e^{i q_k x_n}
const KS_HALF = exp.((KS_Q.^2 .- KS_Q.^4) .* (KS_DT/2))               # e^{Λ dt/2}
const KS_TURN = [cis(2π*k*j/128) for k in 1:KS_K, j in 0:127]          # e^{i q_k σ_j}, σ_j = jL/128
const KS_WORK = (zeros(ComplexF64, KS_K),                              # the nonlinear terms
                 zeros(ComplexF64, KS_K),                              # of the four stages,
                 zeros(ComplexF64, KS_K),
                 zeros(ComplexF64, KS_K),
                 zeros(ComplexF64, KS_K),                              # a stage,
                 zeros(ComplexF64, KS_K),                              # the modes of the state,
                 zeros(KS_N))                                          # and u² on the points

# ---- the nonlinear term -(i q_k/2) (u²)_k of the modes a, into n ----
function ks_nonlinear!(n, a)
    g = KS_WORK[7]
    for p in 1:KS_N
        s = 0.0
        for k in 1:KS_K
            s += 2*real(a[k]*KS_E[k, p])
        end
        g[p] = s^2
    end
    for k in 1:KS_K
        c = zero(ComplexF64)
        for p in 1:KS_N
            c += g[p]*conj(KS_E[k, p])
        end
        n[k] = -0.5im*KS_Q[k]*c/KS_N
    end
    return n
end

# ---- one step of the classical scheme on e^{-Λt} a, in place ----
function ks_step!(a)
    n₁, n₂, n₃, n₄, y = KS_WORK
    E = KS_HALF
    ks_nonlinear!(n₁, a)
    y .= E .* (a .+ (KS_DT/2) .* n₁)
    ks_nonlinear!(n₂, y)
    y .= E .* a .+ (KS_DT/2) .* n₂
    ks_nonlinear!(n₃, y)
    y .= E .* E .* a .+ KS_DT .* E .* n₃
    ks_nonlinear!(n₄, y)
    a .= E .* E .* a .+ (KS_DT/6) .* (E .* E .* n₁ .+ 2 .* E .* (n₂ .+ n₃) .+ n₄)
    return a
end

# ---- one sample ahead, ten steps of 0.025, on the modes of u ----
function ks_advance!(u)
    a = KS_WORK[6]
    for k in 1:KS_K
        a[k] = complex(u[2k - 1], u[2k])
    end
    for _ in 1:10
        ks_step!(a)
    end
    for k in 1:KS_K
        u[2k - 1] = real(a[k])
        u[2k]     = imag(a[k])
    end
    return u
end

# ---- the mode k of the state u, and that of R u ----
mode(u, k)      = complex(u[2k - 1], u[2k])
reflected(u, k) = -conj(mode(u, k))

# ---- the least ‖S(σ) Rʳ b - a‖/‖a‖ over σ and r, with the σ and the r that reach it ----
function ks_distance(a, b)
    na = 2*sum(abs2, a)
    nb = 2*sum(abs2, b)

    # ---- the coefficients of f_r, G_k = 2 conj(a_k) (Rʳ b)_k, for r = 0 and r = 1 ----
    G = (ntuple(k -> 2*conj(mode(a, k))*mode(b, k), Val(KS_K)),
         ntuple(k -> 2*conj(mode(a, k))*reflected(b, k), Val(KS_K)))

    # ---- f_r = Re Σ_k G_k e^{i q_k σ} at 128 shifts, the largest over the shifts and r ----
    f = -Inf
    σ = 0.0
    r = 0
    for s in 0:1, j in 0:127
        g = 0.0
        for k in 1:KS_K
            g += real(G[s + 1][k]*KS_TURN[k, j + 1])
        end
        if g > f
            f = g
            σ = KS_L*j/128
            r = s
        end
    end

    # ---- Newton's method on f_r from the largest, f′ = -Σ q_k Im(G_k e^{i q_k σ}) and
    #      f″ = -Σ q_k² Re(G_k e^{i q_k σ}) ----
    x = σ
    for _ in 1:4
        f′ = 0.0
        f″ = 0.0
        for k in 1:KS_K
            g   = G[r + 1][k]*cis(KS_Q[k]*x)
            f′ -= KS_Q[k]*imag(g)
            f″ -= KS_Q[k]^2*real(g)
        end
        f″ < 0 ||
            break
        x -= f′/f″
    end

    # ---- f_r at the refined shift, kept if it is larger ----
    g = 0.0
    for k in 1:KS_K
        g += real(G[r + 1][k]*cis(KS_Q[k]*x))
    end
    if g > f
        f = g
        σ = mod(x, KS_L)
    end

    return (sqrt(max(na + nb - 2f, 0.0)/na), σ, r)
end

# ---- S(σ) Rʳ b ----
function ks_shift(b, σ, r)
    u = similar(b)
    for k in 1:KS_K
        c         = (r == 0 ? mode(b, k) : reflected(b, k))*cis(KS_Q[k]*σ)
        u[2k - 1] = real(c)
        u[2k]     = imag(c)
    end
    return u
end

# ---- √2 |a_k|, unchanged by the shift and the reflection, of norm ‖a‖ ----
ks_invariant(u) = [sqrt(2)*hypot(u[2k - 1], u[2k]) for k in 1:KS_K]
