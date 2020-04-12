include("examples/example_adjoint_nonuniform.jl");

q₀ = [ 0.60 * x for x in qₙ];
ψ₀ = zeros(N+1);
ψl = zeros(M+1);
ψr = zeros(M+1);
S = 1000;
β = 0.001;
w = 0.000051;                                                               # Априорный параметр апроксимации дельта-функции
q_final, Js, Qs = minimize(q₀, u₀, ulₘ, urₘ, ψ₀, ψl, ψr, Xₙ, N, Tₘ, M, ε, f1_data, f2_data, S = S, β = β, w = w, create_mesh=mshfrm);

frames = [1:20; 21:div(S,50):S];    # Первые двадцать без пропусков
frames = collect(1:div(S, 50):S);   # С пропусками, чтобы всего было 50 кадров
make_minimzation_gif(Js, Qs, qₙ, Xₙ, name = "Minimization_nonuniform.gif", frames_to_write = frames, β = β)
