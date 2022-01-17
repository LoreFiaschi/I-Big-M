include("../src/I_Big_M.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")

using .BAN

# max -x1 + x_2
#
# -x1+ x2 <=  0    y=0
# -x1+ x2 <=  0    y=1
# -x1-2x2 <= -2

function load_param(experiment)

	# standard
	if experiment == 0
		return [-1;  1;  1; -1; zeros(Ban,10)];
	# min l1 norm
	else
		return [-1-η;  1-η;  1-η; -1-η; zeros(Ban,10)];
	end
end

experiment = 1;
M = α;

A = [I -I -I zeros(Ban,4,2);
     zeros(Ban,1,12) ones(Ban,1,2);
	 zeros(Ban,8,4) I M.*[-ones(Ban,4,1) zeros(Ban,4,1); zeros(Ban,4,1) -ones(Ban,4,1)]];

A_dom = [zeros(Ban,3,4) [-1  1  1 -1 zeros(Ban, 1, 4)  0  0;
                         zeros(Ban, 1, 4) -1  1  1 -1  0  0;
                         zeros(Ban, 1, 4) -1  1 -2 -2  0  2]];

A = [A; A_dom];

b = [zeros(Ban,4); 1; zeros(Ban, 11)];

t = [zeros(Int64,5); -ones(Int64,11)];

c = load_param(experiment);

tol = 1e-5;

obj, x, basis, iter = I_Big_M(A, b,c, t, tol=tol, verbose=false, genLatex=false, showprogress=true);

print("\tSolution: ");
println([x[1]-x[2]; x[3]-x[4]]);
print("\tDisjoint flag: "); println(x[end-1:end]);
println("")

nothing
