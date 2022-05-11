include("modify.jl")
#include("NA_Simplex.jl")
include("NA_Simplex_fast.jl")

using LinearAlgebra

function I_Big_M(A::AbstractMatrix{T}, b::AbstractVector{T}, c::AbstractVector{T}, t::Vector{Int};
					tol::Real, verbose::Bool=false, genLatex::Bool=false,
					showprogress::Bool=false) where T <: Number

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
	any(x->x<0, b) && throw(ArgumentError("Entries of vector b must be non-negative"));
	
	# check wether any c[i] is infinite
	any(x->degree(x)>=1, c) && throw(ArgumentError("Entries of vector c must be non-infinite"));

	_A,_b,_c,initial_basis = modify(A,b,c,t);
	#obj, x, basis, iter = na_simplex(_A,_b,_c,initial_basis,tol,verbose,genLatex,showprogress);
	obj, x, basis, iter = na_simplex(_A,_b,_c,initial_basis,tol); # for fast na_simplex

	# fixing of numerical instabilities in infinitesimal components
	xB = _A[:, basis]\_b;
	x[basis] = xB;
	
	return obj, denoise(x[1:size(A,2)], tol), basis, iter;
end
