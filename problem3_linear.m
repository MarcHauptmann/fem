function problem3_linear
% löst -u'' = 1 mit linearen Basisfunktionen

a = 0;
b = 1;
ul = 0;
ur = 0;

n = 129;
h = 1/(n-1);
elements = linspace(a, b, n);

basis = {@(x) 1-x, @(x) x};

% lokale Matrizen
k = [ 1 -1
     -1  1];

% node ordering
B = zeros(n-1, 2);

parfor i=1:n-1
    B(i, :) = [i, i+1];
end

% Assemblierung
tic

A = createMatrix(k, B, 1/h);
fb = createVector(basis, B, elements, @(x) 1);

% Randbedingungen beachten
[A, fb] = dirichletBoundary(A, fb, [1, n], [ul, ur]);

% lösen
[Dinv, L, R] = jacobiDecompose(A);

u = zeros(length(fb), 1);

while norm(A*u - fb) > 1e-4
    u = Dinv*(fb-(L+R)*u);
end

toc

% exakte Lösung
x = linspace(a, b, 2000);
y = -1/2*(x-1).*x;

% Basisfunktionen
Phi = createBasis(basis, elements, x, B);

% plotten
plot(x, Phi*u, x, y);
legend('Näherung', 'exakt', 'location', 'Best');

fprintf('Fehler: %s\n', norm(Phi*u - y'));