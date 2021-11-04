include("../src/I_Big_M.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")

using .BAN

function load_params(experiment)
	
	if experiment==1
		m = 3;  # number of possible locations for warehouses
		n = 5;  # number of customers
		M = α;

		C = [ 500, 500, 500];  # capacity of warehouses
		f = [ 1000, 1000, 1000 ];  # fixed costs of warehouses

		c = [ 4 5 6 8 10;    # fixed cost associated with each shipping from warehouse (row) to customer (column)
			  6 4 3 5 8;
			  9 7 4 3 4];

		D = [ 80, 270, 250, 160, 180];

		return m, n, M, C, f, c, D, 2*m;

	elseif experiment==2
		m = 2;  # number of possible locations for warehouses
		n = 2;  # number of customers
		M = α;

		C = [ 10, 10];  # capacity of warehouses
		f = [ 100, 1 ];  # fixed costs of warehouses

		c = [ 4 4;  # fixed cost associated with each shipping from warehouse (row) to customer (column)
			  4 4];

		D = [ 5, 5 ];  # demand for each customer

		return m, n, M, C, f, c, D, 2*m;
	end

	return;
end

function initialize_problem(experiment)
	
	m, n, M, C, f, c, D, regions = load_params(experiment);

	if length(C) != m || length(D) != n || length(f) != m || size(c) != (m,n)

		throw(ArgumentError("Input data are incorrect"));
		return;
	end

	A = zeros(4*m*n+2*m+n, 3*m*n+2*m);
	A = convert(Matrix{Ban}, A);

	j = m*n +1;

	for i in 1:m*n   # xij = xij1 + xij2

		A[i, i] = 1;
		A[i, j] = -1;
		A[i, j+1] = -1;

		j +=2;

	end

	row_index = m*n+1;
	bx = zeros(Ban, m*n);
	tx = zeros(Int64, m*n);

	j = m*n+1;
	h = 3*m*n+1;

	for (i, counter) in zip(row_index:row_index+m*n-1, 1:m*n)  # xij1 - yi1M <= 0

		A[i, j] = 1;
		A[i, h] = -M;

		j += 2;

		if counter % n == 0

		      h += 2;
		end
	end

	row_index +=  m*n;
	bx1m = zeros(Ban, m*n);
	tx1m = -ones(Int64, m*n);

	j = m*n +1;
	h = 3*m*n+1;

	for (i, counter) in zip(row_index:row_index+m-1,1:m)    # sumj xij1 - Ci*yi1 <= 0

		for col in 1:n

		      A[i, j] = 1;
		      j += 2;

		end

		A[i, h] = -C[counter];

		h += 2;
	end

	row_index += m;
	bc = zeros(Ban, m);
	tc = -ones(Int64, m);


	j = 1;

	for (i, counter) in zip(row_index:row_index+n-1, 1:n)    # sumi xij == Dj

		for col in 1:m

		      A[i, j] = 1;
		      j += n;
		end

		j = 1 + counter;
	end

	row_index += n;
	bd = D;
	td = zeros(Int64, n);

	j = m*n +2;

	for i in row_index:row_index+m*n-1  # xij2 = 0

		A[i, j] = 1;
		j += 2;
	end

	row_index += m*n;
	b2 = zeros(Ban, m*n);
	t2 = zeros(Int64, m*n);

	j = 3*m*n+1;

	for i in row_index:row_index+m-1  # yi1 + yi2 = 1

		A[i, j] = 1;
		A[i, j+1] = 1;
		j += 2;
	end

	row_index += m;
	by = ones(Ban, m);
	ty = zeros(Int64, m);


	j = m*n+1;
	h = 3*m*n+1;
	l = 1;

	for (i, counter) in zip(row_index:row_index+m*n-1, 1:m*n)  # xij1 - yi1 dj <= 0

		A[i, j] = 1;
		A[i, h] = -D[l];

		j += 2;
		l += 1;

		if counter % n == 0

		      h += 2;
		      l = 1;
		end
	end

	bb = zeros(Ban, m*n);
	tb = -ones(Int64, m*n);


	cost_x = reshape(c', m*n);
	cost_y = [ i % 2 == 1 ? f[floor(Int, i/2)+1] : 0 for i in 1:m*2];

	#b =  ones(Ban, size(A, 1), 1) .* hcat(bx, bx1m, bc, bd, b2, by, bb)';
	b = [bx; bx1m; bc; bd; b2; by; bb]
	cost_vector = [-cost_x; zeros(Ban, 2*m*n); -cost_y];
#	cost_vector = ones(Ban, size(A, 2), 1) .* hcat( -cost_x, zeros(1, 2*m*n), -cost_y)';

	t = [tx; tx1m; tc; td; t2; ty; tb];

	return A, b, cost_vector, t, regions;
end

experiment = 2;

A, b, c, t, regions = initialize_problem(experiment);
tol =1e-5;
obj, x, basis, iter = I_Big_M(A, b,c, t, tol=tol, verbose=false, showprogress=true);

print("\tSolution: "); println(x);
print("\tDisjoint flag: "); println(x[end-regions+1:end]);
println("")


