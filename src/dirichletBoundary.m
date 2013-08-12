function [ Anew, bnew] = dirichletBoundary(A, b, ul, ur, leftIndex, rightIndex)
%NEUMANNBOUNDARY erzeugt ein System mit Neumann-Randbedingung
    n = length(b);

    if nargin <= 4
        leftIndex = 1;
        rightIndex = n;
    end
    
    % rechte Seite
    bnew = b;
    
    bnew = bnew - A(:,leftIndex)*ul;
    bnew(leftIndex) = ul;

    bnew = bnew - A(:,rightIndex)*ur;
    bnew(rightIndex) = ur;
    
    % Matrix
    Anew = A;
    
    for index = [leftIndex, rightIndex]
        Anew(index,:) = zeros(1,n);
        Anew(:,index) = zeros(n,1);
        Anew(index, index) = 1;
    end
end

