function fem
a = 0;
b = 3;
ul = 0;
ur = 0;
n = 31;
x = linspace(a,b,n);
h = (b-a)/(n-1)

f = @(x) x.^2;

phi(1) = LinearBasisFunction(-1,0,1);
phi(2) = LinearBasisFunction(0,1,2);

for i=1:n
    phiGlob(i) = LinearBasisFunction(x(i)-h,x(i),x(i)+h);
end

l = length(phi);
m = zeros(l);
k = zeros(l);

for i=1:l    
    for j=1:l
        product = @(x) phi(i).evaluate(x) .* phi(j).evaluate(x);
        diffProduct = @(x) phi(i).evaluateDerivate(x) .* phi(j).evaluateDerivate(x);
        
        m(i,j) = quadgk(product, 0, 1)*h;
        k(i,j) = quadgk(diffProduct, 0, 1)/h;
    end
end

bvector = zeros(n,1);
M = zeros(n,n);
K = zeros(n,n);

for i=1:n    
    bvector(i) = quadgk(@(x) f(x) .* phiGlob(i).evaluate(x), a, b);
        
    for j = 1:n
        if i-j == 1
            M(i,j) = m(1,2);
            K(i,j) = k(1,2);
        elseif i-j == 0
            if i == 1 || i == n
                M(i,j) = m(1,1); 
                K(i,j) = k(1,1);
            else
                M(i,j) = m(1,1) + m(2,2);
                K(i,j) = k(1,1) + k(2,2);
            end
        elseif i-j == -1
            M(i,j) = m(1,2);
            K(i,j) = k(1,2);
        end
    end
end

K
M
bvector
A = K+M

u(2:n-1)=A(2:n-1,2:n-1)\bvector(2:n-1);
u(1) = 0;
u(n) = 0;

f = -1/(1+exp(3));
y = f*exp(x)+(-f-1)*exp(-x)+1;

plot(x, u, x, y);

end