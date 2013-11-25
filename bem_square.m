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
f = @(x,y) sin(2*x*pi);

% Randwerte
g{1} = @(x, y) 0.1*x;
g{2} = @(x, y) 0.1*ones(size(x));
g{3} = @(x, y) 0.1*x;
g{4} = @(x, y) zeros(size(x));

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

for i=1:n
    width = 1/n;
    offset = width/2;
    
    t = offset+(i-1)*width;
    
    %              x y Seite   start     end
    xi(i,:)     = [t 0   1   t-offset t+offset];
    xi(n+i,:)   = [1 t   2   t-offset t+offset];
    xi(2*n+i,:) = [t 1   3   t-offset t+offset];
    xi(3*n+i,:) = [0 t   4   t-offset t+offset];
end

% Randwerte
for i=1:4*n
    side = xi(i, 3);
    
    u(i) = g{side}(xi(i,1), xi(i,2));
end

% Schleife über die Randpunkte
for i=1:4*n
    xi_x = xi(i, 1);
    xi_y = xi(i, 2);
    
    % Schleife über die Basisfunktionen
    for j = 1:4*n
        side = xi(j, 3);
        a = xi(j, 4);
        b = xi(j, 5);
        
        switch side
            case 1
                H(i,j) = quadgk(@(x) dustar(x, 0, xi_x, xi_y, 0, -1), a, b);
                G(i,j) = quadgk(@(x) ustar(x, 0, xi_x, xi_y), a, b);
            case 2
                H(i,j) = quadgk(@(y) dustar(1, y, xi_x, xi_y, 1, 0), a, b);
                G(i,j) = quadgk(@(y) ustar(1, y, xi_x, xi_y), a, b);
            case 3
                H(i,j) = quadgk(@(x) dustar(x, 1, xi_x, xi_y, 0, 1), a, b);
                G(i,j) = quadgk(@(x) ustar(x, 1, xi_x, xi_y), a, b);
            case 4
                H(i,j) = quadgk(@(y) dustar(0, y, xi_x, xi_y, -1, 0), a, b);
                G(i,j) = quadgk(@(y) ustar(0, y, xi_x, xi_y), a, b);
        end
    end
end

b = zeros(4*n,1);

for i=1:4*n
    xi_x = xi(i, 1);
    xi_y = xi(i, 2);
    
    b(i) = quad2d(@(x,y) f(x,y).*ustar(x, y, xi_x, xi_y), 0, 1, 0, 1);
end

b

I = eye(4*n);
right = (H + 0.5*I) * u - b;

q=G\right;

%% Lösung für innere Punkte
np = 20;

[lx, ly] = meshgrid(linspace(0,1,np), linspace(0,1,np));

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
            
            switch side
                case 1
                    g(k) = quadgk(@(x) ustar(x, 0, xi_x, xi_y), a, b);
                    h(k) = quadgk(@(x) dustar(x, 0, xi_x, xi_y,  0, -1), a, b);
                case 2
                    g(k) = quadgk(@(y) ustar(1, y, xi_x, xi_y), a, b);
                    h(k) = quadgk(@(y) dustar(1, y, xi_x, xi_y,  1,  0), a, b);
                case 3
                    g(k) = quadgk(@(x) ustar(x, 1, xi_x, xi_y), a, b);
                    h(k) = quadgk(@(x) dustar(x, 1, xi_x, xi_y,  0,  1), a, b);
                case 4
                    g(k) = quadgk(@(y) ustar(0, y, xi_x, xi_y), a, b);
                    h(k) = quadgk(@(y) dustar(0, y, xi_x, xi_y, -1,  0), a, b);
            end
        end
        
        z(i,j) = g*q-h*u + quad2d(@(x,y) f(x,y).*ustar(x, y, xi_x, xi_y), 0, 1, 0, 1);
    end
end

toc

surf(lx, ly, z);

end

