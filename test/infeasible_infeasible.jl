include("../src/I_Big_M.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")

using LinearAlgebra
using .BAN

M = α;

A = [I -I -I zeros(Ban,2,2); 				# =
	 zeros(Ban,4,2) I -M.*[1 0; 1 0; 0 1; 0 1];	# <
	 zeros(Ban,1,6) ones(Ban,1,2)];			# =

A_dis = [0 0  1 -1  0  0 -1  0; # >
		 0 0  0  0  1  0  0 -4; # >
		 0 0  0  0  1  0  0 -1; # >
		 0 0  1 -1  0  0  1  0; # <
		 0 0  1  0  0  0 -2  0; # <
		 0 0  0  1  0  0 -2  0; # <
		 0 0  0  0  1  0  0 -3; # <
		 0 0  0  0  1  0  0 -5; # <
		 0 0  0  0  0  1  0 -1; # <
];

A = [A; A_dis];

b = [zeros(Ban,6); 1; zeros(Ban,9)];

t = [zeros(Int64,2); -ones(Int64,4); 0; ones(Int64,3); -ones(Int64,6)];

c = [-1; -1; zeros(Ban,6)];

tol = 1e-5;

obj, x, basis, iter = I_Big_M(A, b,c, t, tol=tol, verbose=false, genLatex=false, showprogress=true);  

print("\tSolution: "); 
println([x[1]; x[2]]);

print("\tFull Solution: "); 
println(x);

print("\tDisjoint flag: "); println(x[end-1:end]);
