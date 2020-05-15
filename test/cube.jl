include("../src/I_Big_M.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
using .BAN, LinearAlgebra

A = [Matrix{Float64}(I,3,3) -Matrix{Float64}(I,3,3)];
    
A = convert(Matrix{Ban}, A);
b = ones(Ban, size(A,2), 1) .* [1, 1, 1, 0, 0, 0];
c =  ones(Ban, size(A,1), 1) .* [one(Ban), one(Ban)>>1, one(Ban)>>2];

tol = (one(Ban)>>(length(one(Ban).num)))*1.e-5;

t = [1, 1, 1];

# N.B. swapped b with c
I_Big_M(A, c, -b, t, eps=tol, verbose = true);

return