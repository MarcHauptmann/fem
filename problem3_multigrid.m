function problem3_multigrid
% löst -u'' + u = x^2 mit linearen Basisfunktionen

a = 0;
b = 1;
ul = 0;
ur = 0.1;

n = 128;
elements = linspace(a, b, n+1);

basis = {@(x) 1-x, @(x) x};

% lokale Matrizen
m = 1/6 * [2 1
           1 2];
k = [ 1 -1
     -1  1];
  
elementCounts = round((2*ones(1,log2(n))).^(1:log2(n)));

Matrices = {};
right = {};

for d=elementCounts
    h = 1/d;
    
    % node ordering
    B = zeros(d, 2);

    parfor i=1:d
        B(i, :) = [i, i+1];
    end

    % Assemblierung
    M = createMatrix(m, B, h);
    K = createMatrix(k, B, 1/h);

    A = M+K;

    fb = createVector(basis, B, linspace(a, b, d+1), @(x) x.^2);

    % Randbedingungen beachten
    [A, fb] = dirichletBoundary(A, fb, [1, d+1], [ul, ur]);
   
    index = log2(d);
    Matrices{index} = A
    right{index} = fb;
end

% lösen
function u=solve(A, b, level)
    dim = length(b);
    u = zeros(dim,1);
    [Dinv, L, R] = jacobiDecompose(A{level});

    if level == 1
        for i=1:20
            u = Dinv*(b-(L+R)*u);
        end
    else
        for i=1:10
            u = Dinv*(b-(L+R)*u);
        end

        r = restriction(dim-1)*(b - A{level}*u);
        u = u + prolongation(dim-1) * solve(A, r, level-1);

        for i=1:20
            u = Dinv*(b-(L+R)*u);
        end
    end
end

u = solve(Matrices, right{log2(n)}, log2(n));

% exakte Lösung
x = linspace(a, b, 1000);
y = (exp(1-x).*(exp(1)*ul - ur - 2*exp(1) + 3))/(exp(2)-1) ...
    + (exp(x).*(-ul + exp(1)*ur - 3*exp(1) + 2))/(exp(2)-1) + x.^2 + 2;

% Basisfunktionen
Phi = createBasis(basis, elements, x, B);

% plotten
plot(x, Phi*u, x, y);
legend('Näherung', 'exakt', 'location', 'Best');

end
