include("packages.jl")

function I_Big_M(A::AbstractMatrix{T},b::AbstractMatrix{T},c::AbstractMatrix{T},t::AbstractVector{Int};
						eps::Number=convert(promote_type(T,Float64),1e-5),verbose::Bool=false,genLatex::Bool=false) where T <: Number

	# The problem form
	#
	#   max c^T x
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
	obj, x, base, iter = na_simplex(_A,_b,_c,initial_base,eps,verbose,genLatex);
	
	return obj, x[1:size(A)[2]], base, iter;
end