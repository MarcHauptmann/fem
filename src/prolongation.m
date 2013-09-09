function P = prolongation(d)
%PROLONGATION Erstellt eine Matrix für die Prolongation
    P = zeros(d+1, d/2+1);
    
    for i=1:d+1
        if mod(i, 2) == 1
            P(i, (i-1)/2+1) = 1;
        else
            index = i/2;
            P(i, index:index+1) = [0.5 0.5];
        end
    end
end

