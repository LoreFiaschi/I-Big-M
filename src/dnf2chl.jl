include("../../ArithmeticNonStandarNumbersLibrary/src/BAN_s2_isbits.jl")
#include("../../ArithmeticNonStandarNumbersLibrary/src/BAN_s3_isbits.jl")
#include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
using DelimitedFiles
using LinearAlgebra
using .BAN

function load_constraints(idx)
	A= readdlm("../data/A$(idx).csv", ',')
	#A = convert(Matrix{Float64}, df)
	b = readdlm("../data/b$(idx).csv", ',')
	#b = convert(Matrix{Float64}, df)

	return A, b
end

function dnf2chl(num_pol)

	M = α;
	#M = 100;

	A, b = load_constraints(1);
	Ab = [A b];
	box = [ones(1,size(A,2)) -M];
	convex = [zeros(size(A,2)); 1];

	for i=2:num_pol
		A, b = load_constraints(i)
		sA1, sA2 = size(A)
		Ab = [Ab zeros(size(Ab,1), sA2+1); zeros(sA1, size(Ab,2)) A -b];
		box = [box zeros(size(box,1), sA2+1); zeros(1, size(box,2)) ones(1, sA2) -M];
		convex = [convex; zeros(sA2); 1];
	end

	A_ch = [Ab; box; convex'];
	b_ch = [zeros(Int64, size(A_ch,1)-1); 1];
	t_ch = b_ch.-1;
	b_ch = convert(Vector{Ban}, b_ch);

	return A_ch, b_ch, t_ch
end
