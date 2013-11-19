function laplace_bem2d
%LAPLACE_BEM2D 

g = @(phi) zeros(size(phi));
%g = @(phi) cos(phi)+sin(phi);
f = @(r, phi) r.*sin(phi)+r.*cos(phi);

%     function val = f(r,phi)
%         r
%         
%         val = r.*sin(phi)+r.*cos(phi);
%     end

%f = @(r, phi) cos(r);
%f = @(r, phi) ones(size(r));

% Anzahl Randpunkte
n = 4;

%% Fundamentallösung
    % Fundamentallösung -> -1/(2pi) * ln |x - xi|
    function y = ustar(xr, xphi, xir, xiphi) 
        x_x = xr .* cos(xphi);
        x_y = xr .* sin(xphi);
        
        xi_x = xir .* cos(xiphi);
        xi_y = xir .* sin(xiphi);
        
        dist = sqrt((x_x-xi_x).^2 + (x_y-xi_y).^2);
        
        y = - 1/(2*pi) * log(dist);
    end

    % Ableitung der Fundamentallösung
    function y = dustar(xr, xphi, xir, xiphi, nphi)
        x_x = xr .* cos(xphi);
        x_y = xr .* sin(xphi);
        
        xi_x = xir .* cos(xiphi);
        xi_y = xir .* sin(xiphi);
        
        dist = (x_x-xi_x).^2 + (x_y-xi_y).^2;
        
        n_x = cos(nphi);
        n_y = sin(nphi);
        
        y = -1./(2*pi*dist) .* ((x_x - xi_x).*n_x + (x_y - xi_y).*n_y);
    end

%% Lösen
eps = 1e-15;

    function b = right(xi_phi, g, f, eps)
        alpha = pi;
        b = zeros(size(xi_phi));

        for i = 1:length(b)
            intdustar = quadgk(@(phi) g(phi).*dustar(1, phi, 1, xi_phi(i), phi), xi_phi(i) + eps, xi_phi(i) + 2*pi - eps);
                    
            intf = quad2d(@(r, phi) r.*f(r, phi) .* ustar(r, phi, 1, xi_phi(i)), 0, 1, 0, 2*pi, 'Singular', true);

            b(i,1) = (1 - alpha/(2*pi))*g(xi_phi(i)) + intdustar - intf;
        end
    end


xi = linspace(0, 2*pi, n+1)';
xi = xi(1:n);

b = right(xi, g, f, eps)

A = zeros(n);

dphi = 2*pi/n;

c = {};

for i=1:n
    c{i} = @(phi) max(1 - abs(mod(phi-xi(i), 2*pi))./dphi, 0); % Hutfunktion
end

% i: Index der Basisfunktion
% j: Index des Randpunktes
for i=1:n
    for j = 1:n
        A(j, i) = quadgk(@(phi) ustar(1, phi, 1, xi(j)).*c{i}(phi), xi(j) + eps, xi(j) + 2*pi - eps);
    end
end

A

du = A\b

%% Plotten
nphi = n+1;
%nphi = 100;
nr = 10;

[phi, r] = meshgrid(linspace(0,2*pi,nphi), linspace(0, 1, nr));

x = r.*sin(phi);
y = r.*cos(phi);

alpha = pi;

z = zeros(nr, nphi);

for i = 1:nr
    for j = 1:nphi
        xi_phi = phi(j);
        xi_r = r(i);
        
        intf = quad2d(@(r, phi) r.*f(r, phi) .* ustar(r, phi, xi_r, xi_phi), 0, 1, 0, 2*pi, 'Singular', true);
        
        intdustar = quadgk(@(phi) g(phi).*dustar(1, phi, xi_r, xi_phi, phi), xi_phi + eps, xi_phi + 2*pi - eps);

        intdu = 0;
        
        for k=1:n
            intdu = intdu + du(k)*quadgk(@(phi) ustar(1, phi, xi_r, xi_phi).*c{k}(phi), xi_phi + eps, xi_phi + 2*pi - eps);
        end

        z(i, j) = 1/(1-alpha/(2*pi)) * (intf - intdustar + intdu);
    end
end

%z = ustar(r, phi, 1, phi(2));

surf(x,y,z);

end

