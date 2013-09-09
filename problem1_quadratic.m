function problem1_quadratic
% löst -u'' + u = 1 mit quadratischen Basisfunktionen
% hier gibt es noch einen Fehler

a = 0;
b = 1;
ul = 0;
ur = 0.1;

n = 10;
h = 1/(n-1);
elements = linspace(a, b, n);

basis  = {@(x) 1/2 * x.*(x-1)
          @(x) 1-x.^2
          @(x) 1/2 * x.*(x+1)};
      
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

fb = createVector(basis, B, elements, @(x) x.^2);

% Randbedingungen beachten
[A, fb] = dirichletBoundary(A, fb, [1, 2*n-1], [ul, ur]);

% lösen
u = A\fb

% exakte Lösung
x = linspace(a, b, 1000);
y = (exp(1-x).*(exp(1)*ul - ur - 2*exp(1) + 3))/(exp(2)-1) ...
    + (exp(x).*(-ul + exp(1)*ur - 3*exp(1) + 2))/(exp(2)-1) + x.^2 + 2;

% Basisfunktionen      
Phi = createBasis(basis, elements, x, B, -1, 1);

% plotten
plot(x, Phi*u, x, y);
legend('Näherung', 'exakt', 'location', 'Best');

end