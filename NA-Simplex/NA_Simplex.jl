function NA_Simplex(A::Matrix{T},b::Array{T,2},c::Array{T,2},B::Array{Int64,1},
						eps::Number=convert(promote_type(T,Float64),1e-5),verbose::Bool=true,genLatex::Bool=true) where T <: Number


	n_constraints, n_variables = size(A);
    N = setdiff(1:n_variables, B);
   
    # Assume rank non-deficient initial base matrix
    xB = inv(A[:,B])*b;
    
    any(x->x<0, xB) && error("Unfeasible problem")

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
            print_latex(x/sum(x));
            print(" \$ & \$ ");
            print_latex((c'*x)[1]);
            println(" \$ \\\\");
            println("\t\\hline");
        elseif verbose
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
            
            if genLatex
                println("\\end{tabular}");
                println("\\label{tab:}")
                println("\\end{table}");
                println("");
            elseif verbose
                println("Optimization completed");
                println("Resume:");
                println("\ttotal iterations: $iter");
                print("\tobjective function: "); println(obj);
            end
            
            return obj, x, B;
        end
        
        d = inv_A_B*A[:,N[k]];
        zz = findall(x->x>eps, d);
        
        if isempty(zz)
            obj = convert(T, Inf);
            xB = NaN;
            println("Problem unbounded");
            return obj, xB, B
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