include("packages.jl")

function Big_M(A::AbstractMatrix{T},b::AbstractMatrix{T},c::AbstractMatrix{T},t::AbstractVector{Int}; M::Real=1000,
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

	_A,_b,_c,initial_base = modify(A,b,c,t,M);
    
    init_b = copy(initial_base)
    
	obj, x, base, iter = na_simplex(_A,_b,_c,initial_base,eps,verbose,genLatex);
    
    feasible = ifelse(obj==Inf, obj, all(z->abs(z)<=eps, view(x, view(init_b, count(z->z<0, t)+1:length(init_b)))));
    
    #println(view(init_b, count(z->z<0, t)+1:length(init_b)));
    #println(x[findall(z->abs(z)>eps, x[view(init_b, count(z->z<0, t)+1:length(init_b))])]);
    println("")
    idx_c = findall(z->abs(z)>eps, c)
    idx_x = findall(z->abs(z)>eps, x)
    idx = setdiff(idx_x, idx_c)
    println(idx)
    println(x[idx])
    println("");
    	
	return obj, x[1:size(A)[2]], base, iter, feasible;
end