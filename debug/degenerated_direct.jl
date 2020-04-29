include("../Na-Simplex/NA_Simplex.jl")
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

b = ones(Ban, size(A,1), 1) .* [1, 0.5, 0.5, 0.5, 0.5, 0.75, 0.25];
#b .+= (rand(Ban, length(b),1).>>1).*rand([-1,1], length(b))

c = ones(Ban, 2, 1);

t = [-1, 1, -1, -1, 1, -1, 1];

base = [1,2];

NA_Simplex(A, b, c, t, base);

return