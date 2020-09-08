function modify(A::AbstractMatrix{T},b::AbstractMatrix{T},c::AbstractMatrix{T},t::AbstractVector{Int},M::Real=-1) where T <: Number
	
	# b is assumed with non-negative entries
	# c is assumed to be at most finite
	
	# information gathering
	
	idx_smaller_equal = findall(x->x<0,t);
	idx_equal = findall(x->x==0,t);
	idx_greater_equal = findall(x->x>0,t);
	
	n_constraints, nx = size(A);
	ns = length(idx_smaller_equal);
	ne = length(idx_equal);
	nr = length(idx_greater_equal);
	
	
	# Creation and initialization of the matrix _A (already in the form _A*x <= _b)
	#		 _					   _
	#		|  Ale	I	0	0	0	|
	# _A = 	|  Aeq	0	I	0	0	|
	#		|_ Age  0	0	I  -I  _|
	
	_A = zeros(T, n_constraints, nx+ns+ne+2*nr);
	_A[1:ns, 1:nx] = copy(A[idx_smaller_equal, :]);
	_A[1:ns, nx+1:nx+ns] = Matrix{T}(I,ns,ns); 
	_A[ns+1:ns+ne, 1:nx] = copy(A[idx_equal, :]);
	_A[ns+1:ns+ne, nx+ns+1:nx+ns+ne] = Matrix{T}(I,ne,ne);
	_A[ns+ne+1:ns+ne+nr, 1:nx] = copy(A[idx_greater_equal, :]);
	_A[ns+ne+1:ns+ne+nr, nx+ns+ne+1:nx+ns+ne+nr] = Matrix{T}(I,nr,nr);
	_A[ns+ne+1:ns+ne+nr, nx+ns+ne+nr+1:nx+ns+ne+nr+nr] = -Matrix{T}(I,nr,nr);
	
	
	# Creation and initialization of the vector _b
	#	 	 _				_
	# _b =  |_ ble	beq	 bge_|
	
	_b = Matrix{T}(undef, n_constraints,1);
	_b[1:ns] = copy(b[idx_smaller_equal]);
	_b[ns+1:ns+ne] = copy(b[idx_equal]);
	_b[ns+ne+1:end] = copy(b[idx_greater_equal]);
	
	
	# Creation and initialization of the vector _c
	#	 	 _			   _
	# _c =  |_ c* 0 -1 -1 0_|	*the entries are shrinked in order to be at most infinitesimal of the first order (Still TODO)
	
	_c = zeros(T, nx+ns+ne+2*nr, 1);
	_c[1:nx] = copy(c);
    
    if M < 0 # NA approach
		num = zeros(SIZE);
		num[1]=-1;
        _c[nx+ns+1:nx+ns+ne+nr] .= -α*α*α # Ban(SIZE-1, num); # -α #-ones(T, ne+nr).*α#(α^(degree(maximum(abs.(c)))+length(one(Ban).num)+1));    
    else  # std approach
        _c[nx+ns+1:nx+ns+ne+nr] .= -M;
    end
	
	return _A,_b,_c,[nx+1;nx+2:nx+ns+ne+nr];

end