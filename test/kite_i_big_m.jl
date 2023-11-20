include("../src/I_Big_M.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")

using LinearAlgebra
using .BAN

A = [ 2  1 -3;
      2  3 -2;
      4  3  3;
      0  0  1;
      1  2  1.0]

b = [90, 190, 300, 10, 70.0]

c = [8, 12, 7]

t = [-1, -1, -1, 0, 1]

tol = 1e-5;

A = convert(Matrix{Ban}, A);
b = convert(Vector{Ban}, b);
c = convert(Vector{Ban}, c);

obj, x, B, iter = I_Big_M(A, b, c, t, tol=tol, verbose=true)

println(obj)
println(x)
println(B)
println(iter)
