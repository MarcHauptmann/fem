function problem5
% löst -u'''' = 1 mit linearen Basisfunktionen

a = 0;
b = 1;
ul = 0;
ur = 0;

n = 4;
h = 1/(n-1);
elements = linspace(a, b, n);

basis = {@(x) 1-x, @(x) x};

% lokale Matrizen
m = 1/6 * [2 1
           1 2];
k = [ 1 -1
     -1  1];

% node ordering
B = zeros(n-1, 2);

parfor i=1:n-1
    B(i, :) = [i, i+1];
end

% Assemblierung
tic

M = createMatrix(m, B, h);
K = createMatrix(k, B, 1/h);

A = [M K'; K zeros(n)];

fb = [zeros(n,1); -createVector(basis, B, elements, @(x) -1)];

% Randbedingungen beachten
[A, fb] = dirichletBoundary(A, fb, [n+1, 2*n], [ul, ur])

% lösen
solution = A\fb;

u = solution(n+1:2*n);
p = solution(1:n);

toc

% exakte Lösung
x = linspace(a, b, 2000);
a = -2;
b = 1;
y = 1/24*(x.^4+a*x.^3+b*x.^2);

% Basisfunktionen
Phi = createBasis(basis, elements, x, B);

% plotten
plot(x, Phi*u, x, y);
legend('u', 'exakte Lösung', 'location', 'Best');
