include("../src/I_Big_M.jl")
include("../src/dnf2chl.jl")

num_pol = 100;

A, b, t = dnf2chl(num_pol);

# paying 1 for the ys should not alter the problem
c = ones(Ban, size(A,2));
tol = 1e-5;

fx, x, _, _ = I_Big_M(A, b, c, t, tol=tol, showprogress=true)

println(fx);
println("");
println(x);
