function [ Anew, bnew] = neumannBoundary(A, b, ul, ur)
%NEUMANNBOUNDARY erzeugt ein System mit Neumann-Randbedingung
    n = length(b);

    bnew = b - A(:,1)*ul-A(:,n)*ur;
    bnew(1) = ul;
    bnew(n) = ur;

    Anew = A;
    Anew(1,2:n) = zeros(1,n-1);
    Anew(:,1) = zeros(n,1);
    Anew(:,n) = zeros(n,1);
    Anew(n,2:n) = zeros(1,n-1);
    Anew(1,1) = 1;
    Anew(n,n) = 1;
end

