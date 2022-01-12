include("../src/I_Big_M.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")

using LinearAlgebra
using .BAN

#
# x1 <= 1 && x1 >= 2 && x2 <= 2 && x2 >= 1
#
# -x1+x2 >= 1 && -x1+x2 <= -1
#

M = α;

A = [I -I -I zeros(Ban,4,2);
	zeros(Ban,1,12) ones(Ban,1,2);
	zeros(Ban,8,4) I M.*[-ones(Ban,4,1) zeros(Ban,4,1); zeros(Ban,4,1) -ones(Ban,4,1)]];

A_dom = [zeros(Ban,6,4) [1 -1 zeros(Ban,1,6) -1 0;
						 1 -1 zeros(Ban,1,6) -2 0;
						 0  0 1 -1 zeros(Ban,1,4) -1 0;
						 0  0 1 -1 zeros(Ban,1,4) -2 0;
						zeros(Ban,1,4) 1 -1 -1 1 0 -1;
						zeros(Ban,1,4) 1 -1 -1 1 0  1]];

A = [A; A_dom];

b = [zeros(Ban,4); 1; zeros(Ban, 14)];

t = [zeros(Int64,5); -ones(Int64,8); -1; 1; 1; -1; 1; -1];

c = (-1).*[1; -1;  1; -1; zeros(Ban,10)];

tol = 1e-5;

obj, x, basis, iter = I_Big_M(A, b,c, t, tol=tol, verbose=false, genLatex=true, showprogress=false); 

print("\tSolution: "); 

println([x[1]-x[2]; x[3]-x[4]]);
print("\tDisjoint flag: "); println(x[end-1:end]);
