include("../src/I_Big_M.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")

using .BAN

M = α
experiment = 2;

function load_params(experiment)
	if experiment==1
		return [1; 1; zeros(Ban, size(A,2)-2)];
	end

	return [1; -1; zeros(Ban, size(A,2)-2)];
end

A = [ # x1   x2   x11   x21   x12   x22   x13   x23   x14   x24   y1   y2   y3   y4
        1    0    -1     0    -1     0    -1     0    -1     0    0    0    0    0;    # x1 - x11 - x12 - x13 - x14 = 0
        0    1     0    -1     0    -1     0    -1     0    -1    0    0    0    0;    # x2 - x21 - x22 - x23 - x24 = 0
        0    0     0     0     0     0     0     0     0     0    1    1    1    1;    # y1 + y2 + y3 + y4 = 1
        0    0     1     0     0     0     0     0     0     0   -M    0    0    0;    # x11 - My1 <= 0
        0    0     0     1     0     0     0     0     0     0   -M    0    0    0;    # x21 - My1 <= 0
        0    0     0     0     1     0     0     0     0     0    0   -M    0    0;    # x12 - My2 <= 0
        0    0     0     0     0     1     0     0     0     0    0   -M    0    0;    # x22 - My2 <= 0
        0    0     0     0     0     0     1     0     0     0    0    0   -M    0;    # x13 - My3 <= 0
        0    0     0     0     0     0     0     1     0     0    0    0   -M    0;    # x23 - My3 <= 0
        0    0     0     0     0     0     0     0     1     0    0    0    0   -M;    # x14 - My4 <= 0
        0    0     0     0     0     0     0     0     0     1    0    0    0   -M;    # x24 - My4 <= 0
        0    0     1     1     0     0     0     0     0     0   -1    0    0    0;    # x11 + x21 - y1 >= 0
        0    0     1    -1     0     0     0     0     0     0   -1    0    0    0;    # x11 - x21 - y1 <= 0
        0    0     2    -1     0     0     0     0     0     0   -3    0    0    0;    # 2x11 - x21 - 3y1 <= 0
        0    0    -1     1     0     0     0     0     0     0   -1    0    0    0;    # -x11 + x21 - y1 <= 0
        0    0     1     2     0     0     0     0     0     0   -5    0    0    0;    # x11 + 2x21 - 5y1 <= 0
        0    0     0     0     1    -1     0     0     0     0    0   -5    0    0;    # x12 - x22 - 5y2 >= 0
        0    0     0     0     1    -4     0     0     0     0    0   -5    0    0;    # x12 - 4x22 - 5y2 <= 0
        0    0     0     0     1    -5     0     0     0     0    0    3    0    0;    # x12 - 5x22 + 3y2 <= 0
        0    0     0     0    -1    -3     0     0     0     0    0    14   0    0;    # -x12 - 3x22 + 14y2 >= 0
        0    0     0     0     1    -1     0     0     0     0    0   -7    0    0;    # x12 - x22 - 7y2 <= 0
        0    0     0     0     0     0    -1    -2     0     0    0    0    7    0;    # -x13 - x23 + 7y3 <= 0
        0    0     0     0     0     0    -1     1     0     0    0    0   -2    0;    # -x13 + x23 - 2y3 >= 0
        0    0     0     0     0     0     2    -1     0     0    0    0   3.5   0;    # 2x13 - x23 + 3.5y3 >= 0
        0    0     0     0     0     0     1    -2     0     0    0    0    14   0;    # x13 - 2x23 + 14y3 >= 0
        0    0     0     0     0     0    -1    -1     0     0    0    0    15   0;    # -x13 - x23 + 15y3 >= 0
        0    0     0     0     0     0     0     0    -1    -1    0    0    0   19;    # -x14 - x24 + 19y4 <= 0
        0    0     0     0     0     0     0     0     1    -2    0    0    0    6;    # x14 - 2x24 + 6y4 <= 0
        0    0     0     0     0     0     0     0     1    -1    0    0    0    2;    # x14 - x24 + 2y4 >= 0
        0    0     0     0     0     0     0     0    -1    -2    0    0    0   35;    # -x14 - 2x24 + 35y4 >= 0
        0    0     0     0     0     0     0     0    -1    -1    0    0    0   24;    # -x14 - x24 + 24y4 >= 0

];

b = [zeros(Ban, 2); 1; zeros(Ban, size(A,1)-3)];

c = load_params(experiment);

tol = 1e-5;
t = [zeros(Int64, 3); -ones(Int64, 8); 1; -ones(Int64, 4); 1; -ones(Int64, 2); 1; -ones(Int64, 2); ones(Int64, 4); -ones(Int64, 2); ones(Int64, 3)];

obj, x, basis, iter = I_Big_M(A, b,c, t, tol=tol, verbose=false, genLatex=true, showprogress=false);

print("\tSolution: "); println(x[1:2]);
print("\tDisjoint flag: "); println(x[end-3:end]);
println("")

nothing
