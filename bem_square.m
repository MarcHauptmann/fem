function bem_square
%BEM_SQUARE Lösung der Laplace-Gleichung auf einem Quadrat

% (0,1)     top      (1,1)
%      *------------*
%      |            |
%      |            |
% left |            | right
%      |            |
%      |            |
%      *------------*
% (0,0)    bottom    (1,0)

% rechte Seite
    function int = intF(xi_x, xi_y) 
        f = @(x,y) sin(1*x*pi);
        
        int = quad2d(@(x,y) f(x,y).*ustar(x, y, xi_x, xi_y), 0, 2, 0, 1);
    end

% Punkte des Polygons
border = Polygon([ 0 0
                   2 0
                   2 1
                   0 1 ]);
       
% Normalenvektoren
normal = [  0 -1
            1  0
            0  1
           -1  0 ];
       
% Randwerte
g = { @(x, y) 0.1*x 
      @(x, y) 0.2*ones(size(x)) 
      @(x, y) 0.1*x 
      @(x, y) zeros(size(x)) };

n = 4;

%% Fundamentallösung
    function z = ustar(x, y, xi_x, xi_y)
        dist = sqrt((x-xi_x).^2 + (y-xi_y).^2);
        
        z = -1/(2*pi)*log(dist);
    end

    function z = dustar(x, y, xi_x, xi_y, n_x, n_y)
        diff_x = x-xi_x;
        diff_y = y-xi_y;
        
        dist = diff_x.^2 + diff_y.^2;
        
        z = -(diff_x.*n_x+diff_y.*n_y)./(2*pi*dist);
    end
%% LGS aufstellen und lösen
H = zeros(4*n,4*n);
G = zeros(4*n,4*n);
u = zeros(4*n, 1);
xi = zeros(4*n, 5);

tic

for side = 1:border.edgeCount
    xi_index = (side-1)*n+(1:n);
    
    width = 1/n;
    offset = width/2;
    
    t = offset+((1:n)-1)*width;
    
    xi(xi_index,:) = [border.pos(t, side) side*ones(n,1) t'-offset t'+offset];
end

% Randwerte
for i=1:border.edgeCount*n
    side = xi(i, 3);
    
    u(i) = g{side}(xi(i,1), xi(i,2));
end

% Schleife über die Randpunkte
for i = 1:border.edgeCount*n
    xi_x = xi(i, 1);
    xi_y = xi(i, 2);
    
    % Schleife über die Basisfunktionen
    for j = 1:size(xi, 1)
        side = xi(j, 3);
        a = xi(j, 4);
        b = xi(j, 5);

        n_x = normal(side, 1);
        n_y = normal(side, 2);
        
        H(i,j) = border.speed(side) * quadgk(@(t) dustar(border.x(t, side), border.y(t, side), xi_x, xi_y, n_x, n_y), a, b);
        G(i,j) = border.speed(side) * quadgk(@(t) ustar(border.x(t, side), border.y(t, side), xi_x, xi_y), a, b);
    end
end

b = zeros(border.edgeCount*n,1);

for i=1:4*n
    xi_x = xi(i, 1);
    xi_y = xi(i, 2);
    
    b(i) = intF(xi_x, xi_y);
end

I = eye(border.edgeCount*n);
right = (H + 0.5*I) * u - b;

q=G\right;

%% Lösung für innere Punkte
np = 20;

[lx, ly] = meshgrid(linspace(0,2,np), linspace(0,1,np));

z = zeros(size(lx));

% Randwerte übernehmen
z(1,:) = g{1}(lx(1,:), ly(1,:));
z(np,:) = g{3}(lx(np,:), ly(np,:));
z(:,1) = g{4}(lx(:,1), ly(:,1));
z(:,np) = g{2}(lx(:,np), ly(:,np));
    
% Werte innen berechnen
for i = 2:np-1
    for j = 2:np-1
        xi_x = lx(i, j);
        xi_y = ly(i, j);

        h = zeros(1, 4*n);
        g = zeros(1, 4*n);
        
        for k = 1:4*n
            side = xi(k, 3);
            a = xi(k, 4);
            b = xi(k, 5);
            
            n_x = normal(side, 1);
            n_y = normal(side, 2);
            
            g(k) = quadgk(@(t) border.speed(side) * ustar(border.x(t, side), border.y(t, side), xi_x, xi_y), a, b);
            h(k) = quadgk(@(t) border.speed(side) * dustar(border.x(t, side), border.y(t, side), xi_x, xi_y,  n_x, n_y), a, b);
        end
        
        z(i,j) = g*q-h*u + intF(xi_x, xi_y);
    end
end

toc

surf(lx, ly, z);

end

