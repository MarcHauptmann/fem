function u = multigridSolve(A, b, level)
%MULTIGRIDSOLVE LGS per Multigrid-Verfahren lösen
    dim = length(b);
    
    if level == 1
        u = A\b;
    else
        [Dinv, L, R] = jacobiDecompose(A);
        
        u = zeros(dim,1);
        
        for i=1:100
            u = Dinv*(b-(L+R)*u);
        end

        rest = restriction(dim-2);
        [m,n] = size(rest);
        
        r2 = zeros(m+2, n+2);
        r2(2:m+1,2:n+1) = rest;
        r2(1,1) = 2;
        r2(m+2,n+2) = 2;
        r = 1 * r2;
        
        p = r2';
        p(1:2,1) = [1; 0.5];
        p(n+1:n+2,m+2) = [0.5; 1];
        
        u = u - p * multigridSolve(r*A*p, r*(A*u - b), level-1);

        for i=1:3
            u = Dinv*(b-(L+R)*u);
        end
    end
end

