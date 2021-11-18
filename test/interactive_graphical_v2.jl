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
z5 = xV_*(yB-yC) + yV_*(xC-xB) - xD*(yB-yC) + xC*yB + xB*(yD-yC) -xC*yD;

w1 = yB-yC;
w2 = xC-xB;

M = α;

A_disj_14 = [I zeros(Ban, 2, 8) -I zeros(Ban, 2, 8) [-z1 0; 0 -z2] zeros(Ban, 2, 3); zeros(Ban, 2, 7) I zeros(Ban, 2, 8) -I zeros(Ban, 2, 3) [-z3 0 0; 0 -z4 0]];

A_disj_5 = [zeros(Ban, 1, 4) w1 zeros(Ban, 1, 4) w2 zeros(Ban, 1,4) -w1 zeros(Ban, 1, 4) -w2 zeros(Ban, 1,4) -z5];

A_bounds = [I -M.*[Matrix(I, 5, 5); I; I; I]];

A = [A_disj_14; A_disj_5; A_bounds; zeros(Ban, 1, 20) ones(Ban, 1, 5)];

b = [zeros(Ban, 25); 1];

t = [-1; 1; -1; 1; -1; -ones(Int64, 20); 0];

c = [ones(Ban, 20); zeros(Ban, 5)];

tol = 1e-5;

obj, x, basis, iter = I_Big_M(A, b,c, t, tol=tol, verbose=true, genLatex=false, showprogress=false); 


print("\tSolution: "); 

v1 = xV_ - (sum(x[1:5])-sum(x[11:15]));
v2 = yV_ - (sum(x[6:10])-sum(x[16:20]));

println([v1; v2]);
print("\tDisjoint flag: "); println(x[end-4:end]);
