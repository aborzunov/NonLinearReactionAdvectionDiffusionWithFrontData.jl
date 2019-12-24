using NonLinearReactionAdvectionDiffusionWithFrontData
using LinearAlgebra, ForwardDiff, Printf

u_l(t) = -8#*cos(2*π * t);
u_r(t) = 4# * (1 + sin(2*π * t));
ε = 0.2;
a, b = 0, 1;
t₀, T = 0, 1;
N, M = 40, 80;
h = (b-a)/N;
τ = (T-t₀)/M;
Xₙ = [a  + n*h for n in 0:N];
Tₘ = [t₀ + m*τ for m in 0:M];
q(x) = sin(3 * π * x);
u = zeros(M+1, N+1);

y = u_init.( Xₙ[n] for n in 2:N );
qₙ = [ q(x) for x in Xₙ[2:N] ];
u= solve!(y, Xₙ, Tₘ, N, M, ε, u_l, u_r, qₙ)

make_gif(u, Xₙ, Tₘ; frame_skip = div(M,30), frames_to_write=80);
