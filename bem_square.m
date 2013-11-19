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
u_right = @(x, y) zeros(size(x));
u_top = @(x, y) sin(x*pi);
u_left = @(x, y) zeros(size(x));

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

% Basisfunktionen
for i = 1:4*n
    phi{i} = @(t) t >= xi(i, 4) & t <= xi(i, 5); 
end

% Schleife über die Randpunkte
for i=1:4*n
    xi_x = xi(i, 1);
    xi_y = xi(i, 2);
    
    % Schleife über die Basisfunktionen
    for j = 1:4*n
        side = xi(j, 3);
        
        switch side
            case 1
                H(i,j) = quadgk(@(x) phi{j}(x) .* dustar(x, 0, xi_x, xi_y, 0, -1), 0, 1) ...
                    + (i == j) * (0.5*u_bottom(xi_x, xi_y));
                
                G(i,j) = quadgk(@(x) phi{j}(x) .* ustar(x, 0, xi_x, xi_y), 0, 1);
            case 2
                H(i,j) = quadgk(@(y) phi{j}(y) .* dustar(1, y, xi_x, xi_y, 1, 0), 0, 1) ... 
                    + (i == j) * (0.5*u_right(xi_x, xi_y));
                
                G(i,j) = quadgk(@(y) phi{j}(y) .* ustar(1, y, xi_x, xi_y), 0, 1);
            case 3
                H(i,j) = quadgk(@(x) phi{j}(x) .* dustar(x, 1, xi_x, xi_y, 0, 1), 0, 1) ...
                    + (i == j) * (0.5*u_top(xi_x, xi_y));
                
                G(i,j) = quadgk(@(x) phi{j}(x) .* ustar(x, 1, xi_x, xi_y), 0, 1);
            case 4
                H(i,j) = quadgk(@(y) phi{j}(y) .* dustar(0, y, xi_x, xi_y, -1, 0), 0, 1) ...
                    + (i == j) * (0.5*u_left(xi_x, xi_y));
                
                G(i,j) = quadgk(@(y) phi{j}(y) .* ustar(0, y, xi_x, xi_y), 0, 1);
        end
    end
end

q=G\(H*u);

%% Lösung für innere Punkte
np = 15;

[lx, ly] = meshgrid(linspace(0,1,np), linspace(0,1,np));

z = zeros(size(lx));

for i = 1:np
    z(i,1) = u_bottom(lx(i,1), ly(i,1));
    z(i,np) = u_top(lx(i,np), ly(i,np));
    
    z(1,i) = u_top(lx(1,i), ly(1,i));
    z(np,i) = u_bottom(lx(np,i), ly(np,i));
end
    
for i = 2:np-1
    for j = 2:np-1
        xi_x = lx(i, j);
        xi_y = ly(i, j);

        for k = 1:4*n
            side = xi(k, 3);
            
            switch side
                case 1
                    z(i,j) = z(i,j) + q(k) * quadgk(@(x) phi{k}(x) .* ustar(x, 0, xi_x, xi_y), 0, 1) ...
                            - u(k) * quadgk(@(x) phi{k}(x) .* dustar(x, 0, xi_x, xi_y,  0, -1), 0, 1);
                case 2
                    z(i,j) = z(i,j) + q(k) * quadgk(@(y) phi{k}(y) .* ustar(1, y, xi_x, xi_y), 0, 1) ...
                            - u(k) * quadgk(@(y) phi{k}(y) .* dustar(1, y, xi_x, xi_y,  1,  0), 0, 1);
                case 3
                    z(i,j) = z(i,j) + q(k) * quadgk(@(x) phi{k}(x) .* ustar(x, 1, xi_x, xi_y), 0, 1) ...
                            - u(k) * quadgk(@(x) phi{k}(x) .* dustar(x, 1, xi_x, xi_y,  0,  1), 0, 1);
                case 4
                    z(i,j) = z(i,j) + q(k) * quadgk(@(y) phi{k}(y) .* ustar(0, y, xi_x, xi_y), 0, 1) ...
                            - u(k) * quadgk(@(y) phi{k}(y) .* dustar(0, y, xi_x, xi_y, -1,  0), 0, 1);
            end
        end
    end
end

surf(lx, ly, z);

end

