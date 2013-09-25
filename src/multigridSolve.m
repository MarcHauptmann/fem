function u = multigridSolve(A, b, level)
%MULTIGRIDSOLVE LGS per Multigrid-Verfahren lösen
    dim = length(b);
    
    if level == 1
        u = A\b;
    else
        u = zeros(dim,1);
        [Dinv, L, R] = jacobiDecompose(A);
        
        for i=1:3
            u = Dinv*(b-(L+R)*u);
        end

        rest = restriction(dim-2);
        [m,n] = size(rest);

        r2 = zeros(m+2, n+2);
        r2(2:m+1,2:n+1) = rest;
        r2(1,1) = 2;
        r2(m+2,n+2) = 2;
        r = 0.5 * r2;

        if r(m+1,n+1) == 0
            error('Element 0');
        end
        
        p = r2';
        p(1:2,1) = [1; 0.5];
        p(n+1:n+2,m+2) = [0.5; 1];
            
        for i=1:2
            u = u - p * multigridSolve(r*A*p, r*(A*u - b), level-1);
        end
        
        for i=1:3
            u = Dinv*(b-(L+R)*u);
        end
    end
end

