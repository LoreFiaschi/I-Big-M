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
		sN = denoise(sN, 1e-12)
		
		#print("sN: "); println(sN);
		
		# Gradient Descent
        # k = argmax(sN);
        # k_val = sN[k];
		
		# Bland Rule
		ind_of_pos = findfirst(x->x>eps, sN);
		if ind_of_pos == nothing
			k_val = -1;
		    k = [];
		else
		 	k = ind_of_pos;
			k_val = sN[k];
		end
		
		#print("k_val: "); println(k_val);
		#print("k: "); println(k);
		
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
            
            return obj, x, B, iter;
        end
        
        d = inv_A_B*A[:,N[k]];
        #println(d)
		
		# standard version: findall(x->x>eps, d);
		# the NA counterpart check if any component is > eps and all the previous are at least -eps, i.e., positive
        zz = findall(z->(idx=findfirst(x->x>eps, z.num); idx!=nothing && all(x->x>-eps, z.num[1:idx-1])), d);
        
        if isempty(zz)
            obj = convert(T, Inf);
            x .= NaN;
            verbose && println("Problem unbounded");
            return obj, x, B, iter
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
        x[N[k]] = zero(T);
    end
end