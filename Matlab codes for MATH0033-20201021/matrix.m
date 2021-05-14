function [A,b] = matrix(n,epsi)
% Produces the matrix A and vector b needed for exercise 2 on computational
% exercise sheet 2
A = diag(ones(n,1)) + ...
    epsi.*(diag(ones(n-1,1),-1) + diag(ones(n-1,1),1)) + ...
    (epsi.^2).*(diag(ones(n-2,1),-2) + diag(ones(n-2,1),2));
b = A*ones(n,1);