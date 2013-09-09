function problem2_linear

% Grid
nx = 10;
ny = 10;

points = [-1,  1, 1, -1
          -1, -1, 1,  1 ];
geom(1:2)  = [2, 4];
geom(3:6)  = points(1,:);
geom(7:10) = points(2,:);

dl = decsg(geom');
[p,e,t] = poimesh(dl, nx, ny);

n = size(p, 2);

% node ordering
B = t(1:3,:)';

% Matrix aufstellen
x = 2/nx * [0, 1, 1];
y = 2/ny * [0, 0, 1];
J = (x(2)-x(1))*(y(3)-y(2)) - (x(3)-x(1))*(y(2)-y(1));
a = 1/J * ((x(3)-x(1))^2 + (y(3)-y(1))^2);
b = 1/J * ((y(2)-y(1))*(y(3)-y(1)) + (x(2)-x(1))*(x(3)-x(1)));
c = 1/J * ((x(2)-x(1))^2 + (y(2)-y(1))^2);

S = [ a/2+b+c/2 -a/2-b/2 -b/2-c/2
     -a/2-b/2    a/2      b/2
     -b/2-c/2    b/2      c/2     ];

A = createMatrix(S, B);

% rechte Seite
fb = zeros(n,1);
boundaryPoints = [];

for i=1:n
    topOrBottom = i <= nx+1 | i > (nx+1)^2 - (nx+1);
    leftOrRight = mod(i, nx+1) == 0 | mod(i, nx+1) == 1;
    
    if topOrBottom & leftOrRight
        fb(i) = 1;
        boundaryPoints = [boundaryPoints, i];
    elseif topOrBottom | leftOrRight
        fb(i) = 2;
        boundaryPoints = [boundaryPoints, i];
    else
        fb(i) = 7;
    end
end
fb = 1/6*fb;

% Randwerte
[A, fb] = dirichletBoundary(A, fb, boundaryPoints, zeros(1,length(boundaryPoints)));

% lösen
tic;
u = A\fb;
toc;

% Lösung berechnen
f = {@(x,y) 1-x-y
     @(x,y) x
     @(x,y) y};

% plotten
pdeplot(p, e, t, 'xydata', u, 'zdata', u);

end

