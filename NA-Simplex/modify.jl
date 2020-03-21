function modify(A::Matrix{T},b::Array{T,1},c::Array{T,1},t::Array{Int64,1}) where T = G_Scalar # T <: G_Scalar
	
	# information gathering
	
	idx_smaller_equal = findall(x->x<0,t);
	idx_equal = findall(x->x==0,t);
	idx_greater_equal = findall(x->x>0,t);
	
	n_constraints, nx = size(A);
	ns = length(idx_smaller_equal);
	ne = length(idx_equal);
	nr = length(idx_greater_equal);
	
	type = typeof(A[0,0]);
	
	
	# Creation and initialization of the matrix _A
	#		 _					   _
	#		|  Ale	I	0	0	0	|
	# _A = 	|  Aeq	0	I	0	0	|
	#		|_ Age  0	0	I	-I _|
	
	_A = Matrix{type}(undef, n_constraints, nx+ns+ne+2*nr);
	_A[1:ns, 1:nx] = deepcopy(A[idx_smaller_equal, :]);
	_A[1:ns, nx+1:nx+ns] = Matrix{type}(I,ns,ns); # SI POTRA' RIDEFINIRE?? andare a vedere il codice sorgente in LinearAlgebra
	_A[ns+1:ns+ne, 1:nx] = deepcopy(A[idx_equal, :]);
	_A[ns+1:ns+ne, nx+ns+1:nx+ns+ne] = Matrix{type}(I,ne,ne);
	_A[ns+ne+1:ns+ne+nr, 1:nx] = -1*deepcopy(A[idx_greater_equal, :]);
	_A[ns+ne+1:ns+ne+nr, nx+ns+ne+1:nx+ns+ne+nr] = Matrix{type}(I,nr,nr);
	_A[ns+ne+1:ns+ne+nr, nx+ns+ne+nr+1:nx+ns+ne+nr+nr] = -1*Matrix{type}(I,nr,nr);
	
	
	# Creation and initialization of the vector _b
	#	 	 _				_
	# _b =  |_ ble	beq	-bge_|
	
	_b = deepcopy(b);
	_b[idx_greater_equal] *= -1;
	
	
	# Creation and initialization of the vector _c
	#	 	 _			 _
	# _c =  |_ c* 0 1 1 0_|	*the entries are shrinked in order to be at most infinitesimal of the first order
	
	_c = g_zeros(nx+nr+ne+2*nr); # DA DEFINIRE 
	_c[1:nx] = deepcopy(c);  # RIDEFINIRE OPERATORE ASSEGNAMENTO PER GROSSNUMBERS IN MODO CHE AUTOMATICAMENTE FACCIA LA DEEPCOPY
	scaling_factor = magnitude(max(c))<<1; # DEFINIRE PRINCIPAL COMPONENT, MAGNITUDE (restituisce G-power massima), MAX, << (= *G1)
	_c >>= scaling_factor;
	_c[nx+ns+1:nx+ns+ne+nr] = g_ones(ne+nr); # DA DEFINIRE
	
	
	return _A,_b,_c,[nx+ns+1;nx+ns+2:nx+ns+ne+nr],scaling_factor;

end