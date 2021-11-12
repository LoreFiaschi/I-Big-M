include("../src/I_Big_M.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")

using .BAN

M = α;
experiment = 2;

function load_params(experiment)
	
	if experiment==1
		return [1; 1; zeros(Ban, size(A, 2)-2)];
	end

	return [-1; 1; zeros(Ban, size(A, 2)-2)];

end 

	# x1 x2  x11 x21 x12 x22 y1  y2
A = [ 0   0   0   0   0   0   1   1; # y1 + y2 = 1
	 -1   0   1   0   1   0   0   0; # -x1 + x11 + x12 = 0
	  0  -1   0   1   0   1   0   0; # -x2 + x21 + x22 = 0
	  0   0   1   0   0   0  -M   0; # x11 - M*y1 <= 0
	  0   0   1   0   0   0  -3   0; # x11 - 3*y1 <= 0
	  0   0   0   1   0   0  -M   0; # x21 - M*y1 <= 0
	  0   0   0   1   0   0  -4   0; # x21 - 4*y1 <= 0
	  0   0   0   0   1   0   0  -M; # x12 - M*y2 <= 0
	  0   0   0   0   1   0   0  -5; # x12 - 5*y2 >= 0
	  0   0   0   0   1   0   0  -9; # x12 - 9*y2 <= 0
	  0   0   0   0   0   1   0  -M; # x22 - M*y2 <= 0
	  0   0   0   0   0   1   0  -4; # x22 - 4*y2 >= 0
	  0   0   0   0   0   1   0  -6; # x22 - 6*y2 <= 0
	];

b = [1; zeros(Ban, size(A, 1)-1)];
c = load_params(experiment);

tol = 1e-5;
t = [zeros(Int64, 3); -ones(Int64, 5); 1; -ones(Int64, 2); 1; -1];

obj, x, basis, iter = I_Big_M(A, b,c, t, tol=tol, verbose=false, genLatex=true, showprogress=false);

print("\tSolution: "); println(x[1:2]);
print("\tDisjoint flag: "); println(x[end-1:end]);
println("")

nothing
