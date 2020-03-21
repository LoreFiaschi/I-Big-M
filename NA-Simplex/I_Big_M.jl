function I_Big_M(A::Matrix{T},b::Array{T,1},c::Array{T,1},t::Array{Int64,1},
						eps::Float64=1e-5,verbose::Bool=true,genLatex::Bool=true) where T = G_Scalar # T <: G_Scalar

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

	_A,_b,_c,initial_base,scaling_factor = modify(A,b,c,t);
	obj, x = NA_Simplex(_A,_b,_c,eps,initial_base,verbose,genLatex);
	
	return obj<<scaling_factor, x[1:size(A,2)];
end