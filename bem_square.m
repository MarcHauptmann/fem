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
H = zeros(4,4);
G = zeros(4,4);
u = zeros(4, 1);

% erste Zeile
u(1) = u_bottom(0.5, 1);

% erstes Element
H(1,1) = quadgk(@(x) dustar(x, 0, 0.5, 0, 0, -1), 0, 1) + 0.5*u_bottom(0.5, 0);
G(1,1) = quadgk(@(x) ustar(x, 0, 0.5, 0), 0, 1);

% zweites Element
H(1,2) = quadgk(@(y) dustar(1, y, 0.5, 0, 1, 0), 0, 1);
G(1,2) = quadgk(@(y) ustar(1, y, 0.5, 0), 0, 1);

% drittes Element
H(1,3) = quadgk(@(x) -dustar(x, 1, 0.5, 0, 0, 1), 1, 0);
G(1,3) = quadgk(@(x) -ustar(x, 1, 0.5, 0), 1, 0);

% viertes Element
H(1,4) = quadgk(@(y) -dustar(0, y, 0.5, 0, -1, 0), 1, 0);
G(1,4) = quadgk(@(y) -ustar(0, y, 0.5, 0), 1, 0);

% zweite Zeile
u(2) = u_right(1, 0.5);

% erstes Element
H(2,1) = quadgk(@(x) dustar(x, 0, 1, 0.5, 0, -1), 0, 1);
G(2,1) = quadgk(@(x) ustar(x, 0, 1, 0.5), 0, 1);

% zweites Element
H(2,2) = quadgk(@(y) dustar(1, y, 1, 0.5, 1, 0), 0, 1) + 0.5*u_right(1, 0.5);
G(2,2) = quadgk(@(y) ustar(1, y, 1, 0.5), 0, 1);

% drittes Element
H(2,3) = quadgk(@(x) -dustar(x, 1, 1, 0.5, 0, 1), 1, 0);
G(2,3) = quadgk(@(x) -ustar(x, 1, 1, 0.5), 1, 0);

% viertes Element
H(2,4) = quadgk(@(y) -dustar(0, y, 1, 0.5, -1, 0), 1, 0);
G(2,4) = quadgk(@(y) -ustar(0, y, 1, 0.5), 1, 0);

% dritte Zeile
u(3) = u_top(0.5, 1);

% erstes Element
H(3,1) = quadgk(@(x) dustar(x, 0, 0.5, 1, 0, -1), 0, 1);
G(3,1) = quadgk(@(x) ustar(x, 0, 0.5, 1), 0, 1);

% zweites Element
H(3,2) = quadgk(@(y) dustar(1, y, 0.5, 1, 1, 0), 0, 1);
G(3,2) = quadgk(@(y) ustar(1, y, 0.5, 1), 0, 1);

% drittes Element
H(3,3) = quadgk(@(x) -dustar(x, 1, 0.5, 1, 0, 1), 1, 0) + 0.5*u_right(0.5, 1);
G(3,3) = quadgk(@(x) -ustar(x, 1, 0.5, 1), 1, 0);

% viertes Element
H(3,4) = quadgk(@(y) -dustar(0, y, 0.5, 1, -1, 0), 1, 0);
G(3,4) = quadgk(@(y) -ustar(0, y, 0.5, 1), 1, 0);

% vierte Zeile
u(4) = u_left(0, 0.5);

% erstes Element
H(4,1) = quadgk(@(x) dustar(x, 0, 0, 0.5, 0, -1), 0, 1);
G(4,1) = quadgk(@(x) ustar(x, 0, 0, 0.5), 0, 1);

% zweites Element
H(4,2) = quadgk(@(y) dustar(1, y, 0, 0.5, 1, 0), 0, 1);
G(4,2) = quadgk(@(y) ustar(1, y, 0, 0.5), 0, 1);

% drittes Element
H(4,3) = quadgk(@(x) -dustar(x, 1, 0, 0.5, 0, 1), 1, 0);
G(4,3) = quadgk(@(x) -ustar(x, 1, 0, 0.5), 1, 0);

% viertes Element
H(4,4) = quadgk(@(y) -dustar(0, y, 0, 0.5, -1, 0), 1, 0) + 0.5*u_right(0, 0.5);
G(4,4) = quadgk(@(y) -ustar(0, y, 0, 0.5), 1, 0);

q=G\(H*u);

%% Lösung für innere Punkte
np = 20;

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
        
        z(i,j) = (   q(1) * quadgk(@(x) ustar(x, 0, xi_x, xi_y), 0, 1) ...
                   + q(2) * quadgk(@(y) ustar(1, y, xi_x, xi_y), 0, 1) ...
                   - q(3) * quadgk(@(x) ustar(x, 1, xi_x, xi_y), 1, 0) ...
                   - q(4) * quadgk(@(y) ustar(0, y, xi_x, xi_y), 1, 0) ...
                ) - ( ...
                     u(1) * quadgk(@(x) dustar(x, 0, xi_x, xi_y,  0, -1), 0, 1) ...
                   + u(2) * quadgk(@(y) dustar(1, y, xi_x, xi_y,  1,  0), 0, 1) ...
                   - u(3) * quadgk(@(x) dustar(x, 1, xi_x, xi_y,  0,  1), 1, 0) ...
                   - u(4) * quadgk(@(y) dustar(0, y, xi_x, xi_y, -1,  0), 1, 0));
    end
end

z

surf(lx, ly, z);

end

