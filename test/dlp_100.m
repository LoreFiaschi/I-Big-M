num_pol = 100;
f_opt = zeros(1,num_pol);
x_opt = zeros(2,num_pol);
c = [-1;-1];
A = cell(1,num_pol);
b = cell(1,num_pol);

for i=1:num_pol
    A{1,i} = readmatrix("../data/A"+i+".csv");
    b{1,i} = readmatrix("../data/b"+i+".csv");
end

options = optimoptions('linprog','Algorithm','dual-simplex');

tic
for i=1:num_pol
   [x_opt(:,i), f_opt(i)] = linprog(c, A{1,i}, b{1,i}, [],[], zeros(1,size(A{1,i},2)), [],options);
end

[~, idx] = min(f_opt);
toc
disp(idx)