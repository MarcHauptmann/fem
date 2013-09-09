function [Dinv, L, R] = jacobiDecompose(A)
%JACOBIDECOMPOSE Berechnet eine Zerlegung für das Jacobi-Verfahren

k=length(A);
Dinv = zeros(k);
L = zeros(k);
R = zeros(k);

for i=1:k
    for j=1:k
        if i == j
            Dinv(i,j) = 1/A(i,j);
        elseif i<j
            R(i,j) = A(i,j);
        else
            L(i,j) = A(i,j);
        end
    end
end

end

