function NA_Simplex(A::Matrix{T},b::Array{T,2},c::Array{T,2},B::Array{Int64,1},
						eps::Number=convert(promote_type(T,Float64),1e-5),verbose::Bool=true,genLatex::Bool=true) where T <: Number


	n_constraints, n_variables = size(A);
    N = setdiff(1:n_variables, B);
   
    # Assume rank non-deficient initial base matrix
    xB = inv(A[:,B])*b;
    
    any(x->x<0, xB) && error("Unfeasible problem")

    x = zeros(T, n_variables);
    x[B] = xB;

    iter = 0;
    while true
        iter +=1;
        
        if verbose
            println("Iteration: $iter");
            println(string("\tB: ", B));
            print("\tSolution: ");
            println(x);
            println("");
        end
        
        inv_A_B = inv(A[:,B]);
        
        y = c[B]'*inv_A_B;
        sN = c[N] - A[:,N]'*y';
        
        k = argmax(sN);
        k_val = sN[k];

        if k_val <= eps
            x[B] = xB;
            obj = c'*x;
            println("Optimization completed");
            println("Resume:");
            println("\ttotal iterations: $iter");
            print("\tobjective function: "); println(obj);
            return obj, x, B, y;
        end
        
        d = inv_A_B*A[:,N[k]];
        zz = findall(x->x>eps, d);
        
        if isempty(zz)
            obj = convert(T, Inf);
            x = NaN;
            y = NaN;
            println("Problem unbounded");
            return obj, x, y
        end
        
        quality = xB[zz]./d[zz];
        ii = argmin(quality);
        theta = quality[ii];
        
        l = zz[ii]
        temp = B[l];
        B[l] = N[k];
        N[k] = temp;
        
        xB -= theta*d;
        xB[l] = theta;
     
        x[B] = xB;
        x[N[k]] = zero(Ban);
    end

end