include("../src/I_Big_M.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")

using .BAN

function load_params(experiment)

	if experiment==1
		return 0.3, 1, 3, 3, 1, 1;
	end

	if experiment==2
		return 2, 1, 3, 3, 1, 1;
	end

	if experiment==3
		return 2, 3, 4, 2, 3, 2;
	end

	return 4, 2, 6, 3, 5, 1;
end

M = α
experiment = 3;

E1, E2,	FC1, FC2, VC1, VC2 = load_params(experiment);
PC = 4
PB = 2
B = 10

# x1 represents the quantity of C produced, while x2 represents the cost
A = [#x1  x2   x11  x21  x12  x22  y1   y2   
      1    0   -1    0   -1    0    0    0 ;  # x1 - x11 - x12 = 0
      0    1    0   -1    0   -1    0    0 ;  # x2 - x21 - x22 = 0
      0    0    0    0    0    0    1    1 ;  # y1 + y2 = 1
      0    0    1    0    0    0   -M    0 ;  # x11 - M y1 <= 0
      0    0    0    1    0    0   -M    0 ;  # x21 - M y1 <= 0
      0    0    0    0    1    0    0   -M ;  # x12 - M y2 <= 0
      0    0    0    0    0    1    0   -M ;  # x22 - M y2 <= 0
      0    0    1    0    0    0  -E1*B  0 ;  # x11 - E1 B y1 = 0
      0    0  -VC1   1    0    0  -FC1   0 ;  # -VC1 x11 + x21 - FC1 y1 = 0
      0    0    0    0    1    0    0 -E2*B;  # x12 - E2 B y2 = 0
      0    0    0    0  -VC2   1    0  -FC2;  # -VC2 x12 + x22 - FC2 y2 = 0
];

b = [zeros(Ban, 2); 1; zeros(Ban, size(A,1)-2)];
c = [PC; -1; zeros(Ban, size(A,2)-2)];
t = [zeros(Int64, 3); -ones(Int64, 4); zeros(Int64, 4)];

tol = 1e-5;

obj, x, basis, iter = I_Big_M(A, b,c, t, tol=tol, verbose=false, genLatex=true, showprogress=false);

#print("\tSolution: "); println(x[1:2]);
print("\tSolution: "); println(x);
print("\tDisjoint flag: "); println(x[end-1:end]);
println("")


