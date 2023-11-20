include("../src/NA_Simplex.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")

using LinearAlgebra
using .BAN

A = [ 2  1 -3 1 0 0 0 0 0;
      2  3 -2 0 1 0 0 0 0;
      4  3  3 0 0 1 0 0 0;
      0  0  1 0 0 0 1 0 0;
      0  0 -1 0 0 0 0 1 0;
     -1 -2 -1 0 0 0 0 0 1.0]

b = [90, 190, 300, 10, -10, -70.0]

#c = [8, 12, 7, 0, 0, 0, 0, 0, 0.0]
c = [8+14*η, 12+10*η, 7+2*η, 0, 0, 0, 0, 0, 0.0]

tol = 1e-5;

A = convert(Matrix{Ban}, A);
b = convert(Vector{Ban}, b);
c = convert(Vector{Ban}, c);

B = [2, 3, 4, 5, 6, 8];

obj, x, B, iter = na_simplex(A, b, c, B, tol)

println(obj)
println(x)
println(B)
println(iter)
