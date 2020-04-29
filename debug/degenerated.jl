include("../I_Big_M.jl")
include("../../ArithmeticNonStandarNumbersLibrary/BAN.jl")
using .BAN

A = [ 0      1;
      1      0;
     -1      1;
      1     -1;
      1      1;
     -0.5    1;
      0.5    1];

A = convert(Matrix{Ban}, A);
A -= (rand(Ban, size(A,1), size(A,2)).>>1).*rand([-1,1],size(A,1), size(A,2));

b = ones(Ban, size(A,1), 1) .* [1, 0.5, 0.5, 0.5, 0.5, 0.75, 0.25];
#b .+= (rand(Ban, length(b),1).>>1).*rand([-1,1], length(b))

c = ones(Ban, 2, 1);

t = [-1, 1, -1, -1, 1, -1, 1];

I_Big_M(A, b, c, t, verbose = true);

return