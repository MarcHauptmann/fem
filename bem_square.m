function bem_square
%BEM_SQUARE Lösung der Laplace-Gleichung auf einem Quadrat

% Punkte des Polygons
% border = Polygon([ 0 0
% %                    2 0
%                    2 1
% %                    1 2
%                    0 1 ]);

k = 10;
phi = linspace(0, (1 - 1/k)*2*pi, k);
points = [cos(phi)' sin(phi)'];
border = Polygon(points);

% Randwerte
ug = @(x,y) zeros(size(x));

n = 1;

% rechte Seite
f = @(x,y) 20*exp(-(x.^2+y.^2));

%% Fundamentallösung
    function int = intF(xi_x, xi_y) 
        int = border.integral(@(x,y) f(x, y).*ustar(x, y, xi_x, xi_y));
    end

    function z = ustar(x, y, xi_x, xi_y)
        dist = sqrt((x-xi_x).^2 + (y-xi_y).^2);
        
        z = -1/(2*pi)*log(dist);
    end

    function z = dustar(x, y, xi_x, xi_y, n_x, n_y)
        diff_x = x - xi_x;
        diff_y = y - xi_y;
        
        dist = diff_x.^2 + diff_y.^2;
        
        dotProduct = diff_x.*n_x+diff_y.*n_y;
        
        if min(dotProduct) > 1e-5
            z = -dotProduct./(2*pi*dist);
        else
            z = zeros(size(x));
        end
    end

%% LGS aufstellen und lösen
H = zeros(border.edgeCount*n,border.edgeCount*n);
G = zeros(border.edgeCount*n,border.edgeCount*n);
xi = zeros(border.edgeCount*n, 5);

tic

for side = 1:border.edgeCount
    xi_index = (side-1)*n+(1:n);
    
    width = 1/n;
    offset = width/2;
    
    t = offset+((1:n)-1)*width;
    
    xi(xi_index,:) = [border.pos(t, side) side*ones(n,1) t'-offset t'+offset];
end

% Vektor mit Randwerten
u = ug(xi(:,1), xi(:,2));

% Schleife über die Randpunkte
for i = 1:border.edgeCount*n
    xi_x = xi(i, 1);
    xi_y = xi(i, 2);
    
    % Schleife über die Basisfunktionen
    for j = 1:size(xi, 1)
        side = xi(j, 3);
        a = xi(j, 4);
        b = xi(j, 5);

        normal = border.normal(side);
        n_x = normal(1);
        n_y = normal(2);

        H(i,j) = border.speed(side) * quadgk(@(t) dustar(border.x(t, side), border.y(t, side), xi_x, xi_y, n_x, n_y), a, b);
        G(i,j) = border.speed(side) * quadgk(@(t) ustar(border.x(t, side), border.y(t, side), xi_x, xi_y), a, b);
    end
end

b = zeros(border.edgeCount*n,1);

for i=1:border.edgeCount*n
    xi_x = xi(i, 1);
    xi_y = xi(i, 2);
    
    b(i) = intF(xi_x, xi_y);
end

I = eye(border.edgeCount*n);
right = (H + 0.5*I) * u - b;

q=G\right;

%% Lösung für innere Punkte
[e, tri] = border.triangulation(4);

z = zeros(1, size(e, 1));

for i = 1:length(z)
    xi_x = e(i, 1);
    xi_y = e(i, 2);

    if not(border.contains(xi_x, xi_y))
        h = zeros(1, border.edgeCount*n);
        g = zeros(1, border.edgeCount*n);

        % über alle Basisfunktionen iterieren
        for k = 1:border.edgeCount*n
            side = xi(k, 3);
            a = xi(k, 4);
            b = xi(k, 5);

            normal = border.normal(side);
            n_x = normal(1);
            n_y = normal(2);

            g(k) = quadgk(@(t) border.speed(side) * ustar(border.x(t, side), border.y(t, side), xi_x, xi_y), a, b);
            h(k) = quadgk(@(t) border.speed(side) * dustar(border.x(t, side), border.y(t, side), xi_x, xi_y,  n_x, n_y), a, b);
        end

        z(i) = g*q-h*u + intF(xi_x, xi_y);
    else
        z(i) = ug(xi_x, xi_y);
    end
end

trisurf(tri, e(:,1), e(:,2), z);

end

