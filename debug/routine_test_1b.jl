include("../src/I_Big_M.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
using .BAN
using Debugger

A = Matrix{Ban}(undef, 5, 3);
A[1,1] = Ban(0, [2,0,0,0]); A[1,2] = Ban(0, [1,0,0,0]); A[1,3] = Ban(0, [-3,0,0,0]);
A[2,1] = Ban(0, [2,0,0,0]); A[2,2] = Ban(0, [3,0,0,0]); A[2,3] = Ban(0, [-2,0,0,0]);
A[3,1] = Ban(0, [4,0,0,0]); A[3,2] = Ban(0, [3,0,0,0]); A[3,3] = Ban(0, [3,0,0,0]);
A[4,1] = Ban(0, [0,0,0,0]); A[4,2] = Ban(0, [0,0,0,0]); A[4,3] = Ban(0, [1,0,0,0]);
A[5,1] = Ban(0, [1,0,0,0]); A[5,2] = Ban(0, [2,0,0,0]); A[5,3] = Ban(0, [1,0,0,0]);

b = Array{Ban,2}(undef, 5, 1);
b[1] = Ban(0, [90,0,0,0]);
b[2] = Ban(0, [190,0,0,0]);
b[3] = Ban(0, [300,0,0,0]);
b[4] = Ban(0, [10,0,0,0]);
b[5] = Ban(0, [70,0,0,0]);

c = Array{Ban,2}(undef, 3, 1);
c[1] = Ban(0, [8,14,1,0]);
c[2] = Ban(0, [12,10,1,0]);
c[3] = Ban(0, [7,2,1,0]);

t = [-1, -1, -1, 0, 1];

I_Big_M(A, b, c, t, Ban(0, [0,0,1e-5,0]));

return