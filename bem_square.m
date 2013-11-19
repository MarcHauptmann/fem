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

% Randwerte
u_bottom = @(x, y) sin(x*pi);
u_right = @(x, y) sin(x*pi);
u_top = @(x, y) sin(x*pi);
u_left = @(x, y) sin(x*pi);

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
    
    switch side
        case 1
            u(i) = u_bottom(xi(i,1), xi(i,2));
        case 2
           u(i) = u_right(xi(i,1), xi(i,2));
        case 3
           u(i) = u_top(xi(i,1), xi(i,2));
        case 4
           u(i) = u_left(xi(i,1), xi(i,2));
    end
end

u

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

q=G\(H+0.5*eye(4*n))*u;

%% Lösung für innere Punkte
np = 15;

[lx, ly] = meshgrid(linspace(0,1,np), linspace(0,1,np));

z = zeros(size(lx));

% Randwerte übernehmen
z(1,:) = u_bottom(lx(1,:), ly(1,:));
z(np,:) = u_top(lx(np,:), ly(np,:));
z(:,1) = u_left(lx(:,1), ly(:,1));
z(:,np) = u_right(lx(:,np), ly(:,np));
    
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
        
        z(i,j) = g*q-h*u;
    end
end

toc

surf(lx, ly, z);

end

