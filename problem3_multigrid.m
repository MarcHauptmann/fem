function problem3_multigrid
% löst -u'' = 1 mit linearen Basisfunktionen

a = 0;
b = 1;
ul = 0;
ur = 0;

n = 128;
elements = linspace(a, b, n+1);
h = 1/n;

basis = {@(x) 1-x, @(x) x};

% lokale Matrizen
k = [ 1 -1
     -1  1];

tic 

% node ordering
B = zeros(n, 2);

parfor i=1:n
    B(i, :) = [i, i+1];
end

% Assemblierung
A = createMatrix(k, B, 1/h);
fb = createVector(basis, B, linspace(a, b, n+1), @(x) 1);

% Randbedingungen beachten
[A, fb] = dirichletBoundary(A, fb, [1, n+1], [ul, ur]);

% lösen
u = multigridSolve(A, fb, ceil(log2(n)));

toc

% exakte Lösung
x = linspace(a, b, 1000);
y = -1/2*(x-1).*x;

% Basisfunktionen
Phi = createBasis(basis, elements, x, B);

% plotten
plot(x, Phi*u, x, y);
legend('Näherung', 'exakt', 'location', 'Best');

fprintf('Fehler: %d\n', norm(Phi*u-y'));

end
