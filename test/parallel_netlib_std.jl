using Distributed
@everywhere include("../src/Big_M.jl")
@everywhere using QPSReader
@everywhere using Logging
@everywhere using SparseArrays

# before launching: export JULIA_NUM_THREADS=8

datasets = readdir("netlib");

penalties = [100, 1000, 10000];

#Threads.@threads for file in datasets

    file = datasets[1];

    problem =   with_logger(Logging.NullLogger()) do
                    readqps("netlib/$file");
                end
    
    reverse = 1;
    problem.objsense == :min && (problem.c .*= -1; problem.c0 *= -1; reverse = -1) # reverse the opt direction
    
    #problem.c .*= -1; problem.c0 *= -1; reverse = -1;

    A = spzeros(problem.ncon, problem.nvar);
    for i in eachindex(problem.avals)
        A[problem.arows[i], problem.acols[i]] = problem.avals[i];
    end
    
    A = [A -A; A -A];
    b = [problem.ucon; problem.lcon];
    c = [problem.c; -problem.c];
    t = [-ones(Int64, problem.ncon); ones(Int64, problem.ncon)];
    
    finite_bounds = findall(x->abs(x)!=Inf, b);
    A = A[finite_bounds,:];
    b = b[finite_bounds];
    t = t[finite_bounds];
    
    upper_idx = findall(x->x!=Inf, problem.uvar);
    z = spzeros(length(upper_idx), 2*problem.nvar);
    for i in eachindex(upper_idx)
        z[i, upper_idx[i]] = 1;
        z[i, upper_idx[i]+problem.nvar] = -1
    end
    A = [A; z];
    b = [b; problem.uvar[upper_idx]];
    t = [t; -ones(Int64, length(upper_idx))];
    
    lower_idx = findall(x->x!=-Inf, problem.lvar);
    z = spzeros(length(lower_idx), 2*problem.nvar);
    for i in eachindex(lower_idx)
        z[i, lower_idx[i]] = 1;
        z[i, lower_idx[i]+problem.nvar] = -1
    end
    A = [A; z];
    b = [b; problem.lvar[lower_idx]];
    t = [t; ones(Int64, length(lower_idx))];
    
    for i in eachindex(b)
        b[i] < 0 && (A[i,:] .*= -1; b[i] *= -1; t[i] *= -1;)
    end
    
    tol = 1.e-5;
    
    for M in penalties
        obj, _, _, iter, feasible = Big_M(A, reshape(b, length(b), 1), reshape(c, length(c), 1), t, M=M, eps=tol);
        open("../output/netlib/$(file[1:end-4]).txt", "a+") do out_file
            println(out_file, "$M $(reverse*obj[1]+problem.c0) $iter $feasible");
        end
    end
#end