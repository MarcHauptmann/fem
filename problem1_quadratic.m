function problem1_quadratic
% löst -u'' + u = 1 auf dem Intervall [0, 1] mit den Randwerten 0
% mit quadratischen Basisfunktionen

a = 0;
b = 1;
ul = 0;
ur = 0;

n = 3;
h = 1/(n-1);
x = linspace(a, b, n);

% lokale Matrizen
m = 1/15 * [ 4  2 -1
            2 16  2
           -1  2  4];
k = 1/6 * [ 7 -8  1
           -8 16 -8
            1 -8  7];

% node ordering
B = zeros(n-1, 3);

for i=1:n-1
    num = 1 + 2*(i-1);
    
    B(i, :) = num:num+2;
end

% Assemblierung
M = createMatrix(m, B, h/2);
K = createMatrix(k, B, 2/h);

A = M+K;

fb = h/6 * [1; 4; 2; 4; 1];

% Randbedingungen beachten
[A, fb] = neumannBoundary(A, fb, ul, ur);

% lösen
u = A\fb;

% exakte Lösung
xExact = linspace(a, b, 20);
factor = -1/(1+exp(1));
y = factor*exp(xExact)+(-factor-1)*exp(-xExact)+1;

% Basisfunktionen   
f = {@(x) 1/2 * x.*(x-1)
     @(x) 1-x.^2
     @(x) 1/2 * x.*(x+1)};
             
Phi = createBasis(f, x, xExact, B, -1, 1)

% plotten
plot(xExact, Phi*u, xExact, y);

end