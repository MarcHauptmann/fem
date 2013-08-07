function problem1_linear
% löst -u'' + u = 1 auf dem Intervall [0, 1] mit den Randwerten 0
% mit linearen Basisfunktionen

a = 0;
b = 1;
ul = 0;
ur = 0;

n = 5;
h = 1/(n-1);
x = linspace(a, b, n);

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
M = createMatrix(m, B, h);
K = createMatrix(k, B, 1/h);

A = M+K;

fb = h/2 * [1; 2 * ones(n-2, 1); 1];

% Randbedingungen beachten
[A, fb] = neumannBoundary(A, fb, ul, ur);

% lösen
u = A\fb;

% exakte Lösung
xExact = linspace(a, b, 1000);
factor = -1/(1+exp(1));
y = factor*exp(xExact)+(-factor-1)*exp(-xExact)+1;

% Basisfunktionen
Phi = zeros(length(xExact), n);

for i=1:n
    phi = LinearBasisFunction(x(i)-h, x(i), x(i)+h);
    
    Phi(:,i) = phi.evaluate(xExact);
end

% plotten
plot(xExact, Phi*u, xExact, y);