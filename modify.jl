function modify(A::Matrix{T},b::Array{T,2},c::Array{T,2},t::Array{Int64,1}) where T <: Number
	
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
	_A[ns+ne+1:ns+ne+nr, nx+ns+ne+nr+1:nx+ns+ne+nr+nr] = -1*Matrix{T}(I,nr,nr);
	
	
	# Creation and initialization of the vector _b
	#	 	 _				_
	# _b =  |_ ble	beq	-bge_|
	
	_b = copy(b);
	#_b[idx_greater_equal] *= -1;
	
	
	# Creation and initialization of the vector _c
	#	 	 _			 _
	# _c =  |_ c* 0 -1 -1 0_|	*the entries are shrinked in order to be at most infinitesimal of the first order
	
	_c = zeros(T, nx+ns+ne+2*nr, 1);
	_c[1:nx] = copy(c);
	_c[nx+ns+1:nx+ns+ne+nr] = -ones(T, ne+nr).*(magnitude(maximum(c))<<1);
	
	
	return _A,_b,_c,[nx+1;nx+2:nx+ns+ne+nr];

end