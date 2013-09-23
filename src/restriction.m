function R = restriction(d)
%RESTRICTION Erstellt eine Restriktionsmatrix

    k = (d-1)/2;

    R = zeros(k,d);
    
    a = 0.5;
    
    for i = 1:k
        index = (i-1)*2+1;
        R(i, index:index+2) = [a 1 a];
    end
end

