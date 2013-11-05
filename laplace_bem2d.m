function laplace_bem2d
%LAPLACE_BEM2D 

%% Plotten
n = 20;

[phi, r] = meshgrid(linspace(0,2*pi,n), linspace(0, 1, n));

x = r.*sin(phi);
y = r.*cos(phi);

z = cos(x.^2+y.^2);

surf(x,y,z);

end

