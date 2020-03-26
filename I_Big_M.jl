include("packages.jl")

function I_Big_M(A::Matrix{T},b::Array{T,2},c::Array{T,2},t::Array{Int64,1},
						eps::Number=convert(promote_type(T,Float64),1e-5),verbose::Bool=true,genLatex::Bool=true) where T <: Number

	# The problem form
	#
	#   min c^T x
	#
	#	s.t. Ax <= b for t < 0
	#		 Ax  = b for t = 0
	#		 Ax >= b for t > 0
	#		  x >= 0
	#
	#	Assuming b >= 0
	
	# check whether b >= 0
	any(x->x<0, b) && error("b must be a vector of non-negative entries")

	_A,_b,_c,initial_base = modify(A,b,c,t);
	obj, x = NA_Simplex(_A,_b,_c,initial_base,eps,verbose,genLatex);
	
	return obj, x[1:size(A,2)];
end