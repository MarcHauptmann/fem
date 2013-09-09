function R = restriction(d)
%RESTRICTION Erstellt eine Restriktionsmatrix

    R = zeros(d/2+1,d+1);
    k = d/2+1;
    
    for i=1:k
        if i == 1
            R(i,1:2) = [1 0.5];
        elseif i == k
            R(i,d:d+1) = [0.5 1];
        else
            R(i,2*(i-1):2*(i-1)+2) = [0.5 1 0.5];
        end
    end
end

