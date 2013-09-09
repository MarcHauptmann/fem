function [ Anew, bnew] = dirichletBoundary(A, b, indices, values)
%NEUMANNBOUNDARY erzeugt ein System mit Neumann-Randbedingung
    n = length(b);
    
    % rechte Seite
    bnew = b;
    
    for i=1:length(indices)
        index = indices(i);
        bnew = bnew - A(:,index)*values(i);
        bnew(index) = values(i);
    end
    
    % Matrix
    Anew = A;
    
    for index = indices
        Anew(index,:) = zeros(1,n);
        Anew(:,index) = zeros(n,1);
        Anew(index, index) = 1;
    end
end

