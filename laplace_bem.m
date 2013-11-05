function laplace_bem
%LAPLACE_BEM Löst u''=x^2

%% Vorgaben
% Intervall
a = 0;
b = 2;

% Randwerte
ua = 0;
ub = -1;

% rechte Seite
f = @(x) x.^2;

%% Fundamentallösung
H = @(x) x >= 0;

ustar = @(x, xi) 0.5 * abs(x - xi);
dustar = @(x, xi) -0.5 + H(x - xi);

%% Unbekannte berechnen
eps = 1e-20;

A = [ -ustar(a,a+eps) ustar(b, a) 
      -ustar(a,b) ustar(b, b-eps) ]
  
fb = [ quadgk(@(x) f(x).*ustar(x, a), a, b) + (ub*dustar(b, a) - ua*dustar(a, a+eps)) - ua
       quadgk(@(x) f(x).*ustar(x, b), a, b) + (ub*dustar(b, b-eps) - ua*dustar(a, b)) - ub ]
   
du = A\fb

dua = du(1);
dub = du(2);

%% numerische Lösung
n = 100;

xi = linspace(a, b, n);
y = zeros(1, n);

for i = 2:n-1
    xpos = xi(i);
    
    y(i) = quadgk(@(x) f(x).*ustar(x, xpos), a, b) ...
        + (ub*dustar(b, xpos) - ua*dustar(a, xpos)) ...
        - (dub*ustar(b, xpos) - dua*ustar(a, xpos));
end

y(1) = ua;
y(n) = ub;

%% exakte Lösung
c = [a 1; b 1]\[ua - 1/12*a^4; ub - 1/12*b^4];

yex = 1/12*xi.^4+c(1)*xi+c(2);

%% Plotten
plot(xi, y, xi, yex);
legend('Näherung', 'exakt', 'location', 'best');

end

