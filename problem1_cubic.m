function problem1_cubic
% löst -u'' + u = 1 mit quadratischen Basisfunktionen

a = 0;
b = 1;
ul = 0;
ur = 0.1;

n = 5;
h = 1/(n-1);
elements = linspace(a, b, n);

basis  = {@(x) 1 - 3*x.^2 + 2*x.^3
          @(x) h*(x - 2*x.^2 + x.^3)
          @(x) 3*x.^2 - 2*x.^3
          @(x) h*(-x.^2 + x.^3)};
      
% lokale Matrizen
m = h/420 * [  156   22*h   54   -13*h
               22*h  4*h^2  13*h -3*h^2
               54    13*h   156  -22*h
              -13*h -3*h^2 -22*h  4*h^2 ];  
k = 1/(30*h) * [ 36   3*h   -36   3*h
                 3*h  4*h^2 -3*h -h^2
                -36  -3*h    36  -3*h
                 3*h -h^2   -3*h  4*h^2 ];

% node ordering
B = zeros(n-1, 4);

for i=1:n-1
    num = 1 + 2*(i-1);
    
    B(i, :) = num:num+3;
end

% Assemblierung
M = createMatrix(m, B);
K = createMatrix(k, B);

A = M+K;

fb = createVector(basis, B, elements, @(x) x.^2);

% Randbedingungen beachten
[A, fb] = dirichletBoundary(A, fb, ul, ur, 1, 2*(n-1)+1);

% lösen
u = A\fb;

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

