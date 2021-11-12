using ProgressMeter

function na_simplex(A::AbstractMatrix{T}, b::AbstractVector{T}, c::AbstractVector{T}, B::Vector{Int64},
						tol::Real, verbose::Bool=true, genLatex::Bool=false,
						showprogress::Bool=false) where T <: Number


	# The method implements the revised simplex method in Box 7.1 on page 103 of Chvatal

	# Revised Simplex
	#
	# max  c'*x
	# s.t. Ax = b
	#      x >= 0

	n_constraints, n_variables = size(A);
    N = setdiff(1:n_variables, B);
   
    # Assume rank non-deficient initial base matrix
    xB = (A[:,B]\I)*b; # inversion of A[:,B] by backslash
    
    any(x->x<-tol, xB) && (println(""); true;) && (println(tol); println(""); println(xB[findall(x->x<-tol, xB)]); true;) && error("Unfeasible problem")

    x = zeros(T, n_variables);
    x[B] = xB;
	
	aux_var = map(z->z[1], findall(z->z.p>0, c)); # TODO : careful, it assumes use of BANs
    if genLatex
        println("\\begin{table}[ht]");
        println("\t\\centering");
        println("\t\\caption{}");
        println("\t\\begin{tabular}{|c|c|c|c|}");
        println("\t\\hline");
        println("\t\\textbf{Iter.} & \\textbf{Basis} & \$ \\mathbf{x} \$ & \$ \\tilde{\\mathbf{c}}^T \\mathbf{x} \$ \\\\");
        println("\t\\hline");
    elseif showprogress
		prog = ProgressUnknown("Iteration:");
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
			print_latex(c'*x);
            println(" \$ \\\\");
            println("\t\\hline");
        elseif verbose
            println("Iteration: $iter");
            println(string("\tB: ", B));
            print("\tCost: ")
            println(c'*x)
            println("");
		elseif showprogress
			ProgressMeter.next!(prog);
        end
        
		inv_A_B = A[:,B]\I;
        
        y = c[B]'*inv_A_B;
        sN = c[N] - A[:,N]'*y';
		sN = denoise(sN, tol) # DANGER!! entries can change sign, is it a problem?

		# Bland Rule

		ind_of_pos = findfirst(x->x.num[1]>tol, sN);
		
		if ind_of_pos == nothing
			k_val = -1;
		    k = [];
		else
		 	k = ind_of_pos;
			k_val = sN[k];
		end

        if k_val < 0
            obj = c'*x;
            
            if genLatex
                println("\\end{tabular}");
                println("\\label{tab:}")
                println("\\end{table}");
                println("");
			elseif showprogress
				ProgressMeter.finish!(prog);
            end
			
			print("Optimization completed, ");
			(all(z->z==0, x[aux_var])) ? println("feasible solution found") : println("unfeasible solution found");
			println("Resume:");
			println("\ttotal iterations: $iter");
			print("\tobjective function: "); println(obj);
			println("");
            
            return obj, x, B, iter;
        end
        
        d = denoise(inv_A_B*A[:,N[k]], tol);
		zz = findall(x->x>0, d); # it works thanks to previous denoise
        
        if isempty(zz)
            obj = convert(T, Inf);
            x .= NaN;
            verbose && println("Problem unbounded");
            return obj, x, B, iter
        end
        
        quality = xB[zz]./d[zz];
        ii = argmin(quality); 
        θ = quality[ii];
        
        x[B] -= θ*d 
		x[N[k]] = θ;
		
		l = zz[ii]


        temp = B[l];
        B[l] = N[k];
        N[k] = temp;
		
		x = denoise(x, tol);
		xB = x[B];
    end
end
