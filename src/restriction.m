function R = restriction(d)
%RESTRICTION Erstellt eine Restriktionsmatrix

    k = (d-1)/2;

    R = zeros(k,d);
    
    for i = 1:k
        index = (i-1)*2+1;
        R(i, index:index+2) = [0.5 1 0.5];
    end
end

