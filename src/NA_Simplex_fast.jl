using ProgressMeter

function na_simplex(A::AbstractMatrix{T}, b::AbstractVector{T}, c::AbstractVector{T}, B::Vector{Int64},
						tol::Real) where T <: Number


	# The method implements the revised simplex method in Box 7.1 on page 103 of Chvatal

	# Revised Simplex
	#
	# max  c'*x
	# s.t. Ax = b
	#      x >= 0

	#showprogress = true

	n_constraints, n_variables = size(A);
    N = setdiff(1:n_variables, B);

	# Assume rank non-deficient initial base matrix
	inv_A_B = A[:,B]\I;
    xB = inv_A_B*b; # inversion of A[:,B] by backslash
    
    any(x->x<-tol, xB) && (println(""); true;) && (println(tol); println(""); println(xB[findall(x->x<-tol, xB)]); true;) && error("Unfeasible problem")

    x = zeros(T, n_variables);
    x[B] = xB;
	
	aux_var = map(z->z[1], findall(z->z.p>0, c)); # TODO : careful, it assumes use of BANs
    #if showprogress
		prog = ProgressUnknown("Iteration:");
	#end

	

    iter = 0;
    while true #iter < 44 #
        iter +=1;
        
        #if showprogress
			ProgressMeter.next!(prog);
        #end
        
		#inv_A_B = A[:,B]\I;
        
        y = c[B]'*inv_A_B;
        sN_ = c[N] - A[:,N]'*y';
		sN = denoise(sN_, tol) # DANGER!! entries can change sign, is it a problem?

		# Bland Rule

#		ind_of_pos = findfirst(x->x.num[1]>tol, sN); # General purpose library
		ind_of_pos = findfirst(x->x.num1>tol, sN); # s3 isbits library
		
		if ind_of_pos == nothing
			k_val = -1;
		    k = [];
		else
		 	k = ind_of_pos;
			k_val = sN[k];
		end

        if k_val < 0
            obj = c'*x;
            
			#if showprogress
				ProgressMeter.finish!(prog);
            #end
			
			print("Optimization completed, ");
			(all(z->z==0, x[aux_var])) ? println("feasible solution found") : println("unfeasible solution found");
			println("Resume:");
			println("\ttotal iterations: $iter");
			print("\tobjective function: "); println(obj);
			println("");
            
            return obj, x, B, iter;
        end
        
		d_ = inv_A_B*A[:,N[k]]
        d = denoise(d_, tol);

		zz = findall(x->x>0, d); # it works thanks to previous denoise
        
        if isempty(zz)
            obj = convert(T, Inf);
            x .= NaN;
            verbose && println("Problem unbounded");
            return obj, x, B, iter
        end
        
        quality = xB[zz]./d[zz];
        ii = argmin(quality); 

#        θ = quality[ii];
        
#       x[B] -= θ*d 
#		x[N[k]] = θ;

		l = zz[ii]


        temp = B[l];
        B[l] = N[k];
        N[k] = temp;

		x[N[k]] = 0;

		# reuse of inversion just for computational stability
		inv_A_B = A[:,B]\I;
		x[B] = inv_A_B*b;
		
		x = denoise(x, tol);
		xB = x[B];
    end
end
