function [Dinv, L, R] = jacobiDecompose(A)
%JACOBIDECOMPOSE Berechnet eine Zerlegung für das Jacobi-Verfahren

D = diag(diag(A));

Dinv = diag(1./diag(A));

L = tril(A) - D;
R = triu(A) - D;
end

