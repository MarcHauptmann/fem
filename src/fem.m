function u = fem(x,ul,ur, phi, phiGlob, funcA, funcB, funcC, f)
    n = length(x);

    function m=createLocalMatrix(phi, form)
        m = zeros(length(phi));
        
        for i=1:l
            for j=1:l
                m(i,j) = form(phi(i), phi(j));
            end
        end
    end

    function [A,b] = createGlobalSystem(m,k,n, funcC, f, phiGlob)
        A = zeros(n);
        
        n = length(phiGlob);
        b = zeros(n,1);
        
        for i=1:n
            b(i) = funcC(f, phiGlob(i));
        end
        
        for i=1:n
            for j = 1:n
                pos = i-j;
                
                if pos == 1
                    h = x(i) - x(j);
                        
                    A(i,j) = h * m(2,1) + 1/h * k(2,1);
                elseif pos == 0
                    if i == 1
                        h = x(2) - x(1);
                        
                        A(i,j) = h * m(2,2) + 1/h * k(2,2);
                    elseif i == n
                        h = x(n) - x(n-1);
                        
                        A(i,j) = h * m(1,1) + 1/h * k(1,1);
                    else
                        hl = x(i)-x(i-1);
                        hr = x(i+1)-x(i);
                        
                        A(i,j) =  hl * m(1,1) + hr * m(2,2) + 1/hl * k(1,1) + 1/hr * k(2,2);
                    end
                elseif pos == -1
                    h = x(j) - x(i);
                    
                    A(i,j) = h * m(1,2) + 1/h * k(1,2);
                end
            end
        end
        
        b = b-A(:,1)*ul-A(:,n)*ur;
        b(1) = ul;
        b(n) = ur;
        
        A(1,2:n) = zeros(1,n-1);
        A(:,1) = zeros(n,1);
        A(:,n) = zeros(n,1);
        A(n,2:n) = zeros(1,n-1);
        A(1,1) = 1;
        A(n,n) = 1;
    end

l = length(phi);
m = createLocalMatrix(phi, funcA);
k = createLocalMatrix(phi, funcB);

[A, b] = createGlobalSystem(m,k,n,funcC, f, phiGlob);

u=A\b;

end