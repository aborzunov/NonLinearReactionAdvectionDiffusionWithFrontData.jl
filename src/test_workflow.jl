using NonLinearReactionAdvectionDiffusionWithFrontData
using NonLinearReactionAdvectionDiffusionWithFrontData: DirectProblem
using NonLinearReactionAdvectionDiffusionWithFrontData: u_initial
using NonLinearReactionAdvectionDiffusionWithFrontData: create_mesh
using NonLinearReactionAdvectionDiffusionWithFrontData: solve
using NonLinearReactionAdvectionDiffusionWithFrontData: u_init
using NonLinearReactionAdvectionDiffusionWithFrontData: f
using NonLinearReactionAdvectionDiffusionWithFrontData: solve!
using NonLinearReactionAdvectionDiffusionWithFrontData: make_gif

l(t) = -1#* cos(π*t/2) - 1; # Граничные условия на левой  границе
r(t) = 1 #* sin(π*t/2) + 1;   # Граничные условия на правой границе
a, b = 0, 1
t₀, T = 0, 1
N, M = 1000, 1000
Xn, h = create_mesh(a, b, N)
Tm, tau = create_mesh(t, T, M)
ε = 0.2;
q(x) = 0.5
p = DirectProblem(a, b, t, T, ε, N, M, l, r, q, Xn, Tm, h, tau )

plot( p.ul.(p.Tₘ), label="left")
plot!( p.ur.(p.Tₘ), label="right")
plot!( p.q.(p.Xₙ), label="q")
plot!( u, label="Initial")
plot!( rp(u, 0, p), label="rp")




u = u_initial(p)
answer, history = solve(p)
a = Animation()
for t in 1:5
    plot(p.Xₙ,  history[:,t], ylims=(-2,2))
    frame(a)
end
gif(a,"solution.gif")

jacobian(u,0,p)

out = similar(u);
for n in 1:N-1
    out[n] = (p.ε*π^2*(1 - 2t) + 2t + π * (1 - 2t)^2 * cos(π * p.Xₙ[n+1]) - p.q(p.Xₙ[n+1]) * (1 - 2t) ) * sin(π * p.Xₙ[n+1])
end
plot(out)
plot!(rp(u,0,p))

u_l(t) = -8;
u_r(t) = 4;
ε = 0.2
a, b = 0, 1;
t₀, T = 0, 1;
N, M = 1000, 1000;
h = (b-a)/N;
τ = (T-t₀)/M;
Xₙ = [a  + n*h for n in 0:N];
Tₘ = [t₀ + m*τ for m in 0:M];
q(x) = sin(3 * π * x)
u = zeros(M+1, N+1);

y = u_init.( Xₙ[n] for n in 2:N )
u = solve!(u, y, Xₙ, Tₘ, N, M, ε, u_l, u_r, q)

make_gif(u, Xₙ; frame_skip = 2);
