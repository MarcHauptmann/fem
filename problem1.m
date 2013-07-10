function problem1
a = 0;
b = 3;
ul = 0;
ur = 0;
x = [linspace(0, 1.5, 7), linspace(1.6, 3, 4)]
n = length(x);
phi(1) = LinearBasisFunction(-1,0,1);
phi(2) = LinearBasisFunction(0,1,2);

for i=1:n
    if i == 1
        hl = 1;
    else
        hl = x(i)-x(i-1);
    end
    
    if i == n
        hr = 1;
    else
        hr = x(i+1) - x(i);
    end
    
    phiGlob(i) = LinearBasisFunction(x(i)-hl,x(i),x(i)+hr);
end

    function res = funcM(f1, f2)
        product = @(x) f1.evaluate(x) .* f2.evaluate(x);
        
        res = quadgk(product, 0, 1);
    end


    function res = funcK(f1, f2)
        product = @(x) f1.evaluateDerivate(x) .* f2.evaluateDerivate(x);
        
        res = quadgk(product, 0, 1);
    end

    function res = funcF(f1, f2)
        fun = @(x) f1.evaluate(x) .* f2.evaluate(x);
        
        res = quadgk(fun, a, b);
    end


f = GenericFunction(@(x) ones(1,length(x)));

u = fem(x, ul, ur, phi, phiGlob, @funcM, @funcK, @funcF, f);

factor = -1/(1+exp(3));
y = factor*exp(x)+(-factor-1)*exp(-x)+1;

plot(x, u, x, y);
end