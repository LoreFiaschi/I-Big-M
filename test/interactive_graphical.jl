include("../src/I_Big_M.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")

using LinearAlgebra
using .BAN

xA = 4;
xB = 6;
xC = 4;
xD = 1;
xE = 3;

yA = 1.5;
yB = 2;
yC = 2.5;
yD = 1;
yE = 2;

xV_ = -1.25;
yV_ = 0.25;

z1 = xA-xE+xV_;
z2 = xB-xD+xV_;
z3 = yA-yE+yV_;
z4 = yC-yD+yV_;

M = α;

A_disj = [I zeros(Ban, 2, 6) -I zeros(Ban, 2, 6) [-z1 0; 0 -z2] zeros(Ban, 2, 2); zeros(Ban, 2, 6) I zeros(Ban, 2, 6) -I zeros(Ban, 2, 2) [-z3 0; 0 -z4]];

A_bounds = [I -M.*[Matrix(I, 4, 4); I; I; I]];

A = [A_disj; A_bounds; zeros(Ban, 1, 16) ones(Ban, 1, 4)];

b = [zeros(Ban, 20); 1];

t = [-1; 1; -1; 1; -ones(Int64, 16); 0];

c = [ones(Ban, 16); zeros(Ban, 4)];

tol = 1e-5;

obj, x, basis, iter = I_Big_M(A, b,c, t, tol=tol, verbose=false, genLatex=true, showprogress=false); 


print("\tSolution: "); 

v1 = xV_ - (sum(x[1:4])-sum(x[9:12]));
v2 = yV_ - (sum(x[5:8])-sum(x[13:16]));

println([v1; v2]);
print("\tDisjoint flag: "); println(x[end-3:end]);
