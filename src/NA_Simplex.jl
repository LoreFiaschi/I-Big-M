function na_simplex(A::Matrix{T},b::Array{T,2},c::Array{T,2},B::Array{Int64,1},
						eps::Number=convert(promote_type(T,Float64),1e-5),verbose::Bool=true,genLatex::Bool=false) where T <: Number


	# The method implements the revised simplex method in Box 7.1 on page 103 of Chvatal

	# Revised Simplex
	#
	# max  c'*x
	# s.t. Ax = b
	#      x >= 0

	n_constraints, n_variables = size(A);
    N = setdiff(1:n_variables, B);
   
    # Assume rank non-deficient initial base matrix
    xB = inv(A[:,B])*b;
    
    any(x->x<-eps, xB) && (println(""); true;) && (println(eps); println(""); println(xB[findall(x->x<-eps, xB)]); true;) && error("Unfeasible problem")
    #any(x->x<0, xB) && (println(""); true;) && (println(xB); true;) && (println(b); true;) && error("Unfeasible problem")

    x = zeros(T, n_variables);
    x[B] = xB;
	
	aux_var = map(z->z[1], findall(z->z.p>0, c));
	
	#println(not_aux_var);
	#error()
    
    if genLatex
        println("\\begin{table}[!ht]");
        println("\t\\centering");
        println("\t\\caption{}");
        println("\t\\begin{tabular}{|c|c|c|c|}");
        println("\t\\hline");
        println("\t\\textbf{Iter.} & \\textbf{Base} & \$ \\mathbf{x}^* \$ & \$ \\tilde{\\mathbf{c}}^T \\mathbf{x}^* \$ \\\\");
        println("\t\\hline");
    end

    iter = 0;
    while true	
        iter +=1;
        
        if genLatex
            print("\t$iter & \$ \\{$(B[1])");
            for elem in B[2:end]
            print(",\\, $elem");
            end
            print("\\} \$ & \$ ");
			print_latex(x);
			print(" \$ & \$ ");
			print_latex((c'*x)[1]);
            println(" \$ \\\\");
            println("\t\\hline");
        elseif verbose
            println("Iteration: $iter");
            println(string("\tB: ", B));
            print("\tCost: ")
            println((c'*x)[1])
            #print("\tSolution: ");
            #println(x);
            println("");
        end
        
		#=
		for r in eachrow(A[:,B])
			println(r);
		end
		println("");
		=#
		
        inv_A_B = inv(A[:,B]);
        
        y = c[B]'*inv_A_B;
        sN = c[N] - A[:,N]'*y';
		#print("sN: "); println(sN);
		#println("");
		sN = denoise(sN, eps[1]) # DANGER!! entries can change sign, is it a problem?
		
		#print("sN: "); println(sN);
		#println("");
		
		#not_aux_var_N = findall(z->(z in not_aux_var), N);
		
		# Gradient Descent
        # k = argmax(sN);
        # k_val = sN[k];
		
		# Bland Rule
		#ind_of_pos = findfirst(x->x>eps, sN[not_aux_var_N]); # TODO modify in the external version of >  if considered necessary
		ind_of_pos = findfirst(x->x>eps, sN); # TODO modify in the external version of >  if considered necessary
		if ind_of_pos == nothing
			k_val = -1;
		    k = [];
		else
		 	#k = not_aux_var_N[ind_of_pos];
		 	k = ind_of_pos;
			k_val = sN[k];
		end
		
		#=
		print("k_val: "); println(k_val);
		print("k: "); println(k);
		println("");
		#println("");
		=#
		
		#degree_improvement = findfirst(x->x>0, k_val.num);

        if k_val < 0 #degree_improvement==nothing || any(x->x<=0, (k_val-eps).num[1:degree_improvement])
            #x[B] = xB;
            obj = c'*x;
            
            if genLatex
                println("\\end{tabular}");
                println("\\label{tab:}")
                println("\\end{table}");
                println("");
            end
			
			println("Optimization completed");
			println("Resume:");
			println("\ttotal iterations: $iter");
			print("\tobjective function: "); println(obj);
			println("");
			print("\tauxiliary variables: "); println(denoise(x[aux_var], eps[1]));
			println("");
            
            return obj, x, B, iter;
        end
        
        d = inv_A_B*A[:,N[k]];
        #print("d_pre: "); println(d); println("");
		
		# standard version
		#zz = findall(x->x>eps, d); # TODO transform in the external version of > if considered necessary
		# the NA counterpart check if any component is > eps and all the previous are at least -eps, i.e., positive
        zz = findall(z->(idx=findfirst(x->x>eps[1], z.num); idx!=nothing && all(x->x>-eps[1], z.num[1:idx-1])), d); # Could need of a denoise later
        
        if isempty(zz)
            obj = convert(T, Inf);
            x .= NaN;
            verbose && println("Problem unbounded");
            return obj, x, B, iter
        end
        
        quality = xB[zz]./d[zz];
        ii = argmin(quality); 
        theta = quality[ii];
		
		#=
		print("xB: "); println(xB[zz]);
		print("d: "); println(d[zz]);
		println("");
		print("quality: "); println(quality);
		print("theta: "); println(theta);
		#print("ii: "); println(ii);
		println("");
		println("");
        =#
		
        
        
        x[B] -= theta*d 
		x[N[k]] = theta;
		
		l = zz[ii]
        temp = B[l];
        B[l] = N[k];
        N[k] = temp;
		
		xB = x[B];
		
		#neg = findall(z->z<0, x);
		#neg!=[] && (print("negative entries: "); println(neg); print("value: "); println(x[neg]); println("");)
		#x = denoise(x, eps[1]);
    end
end