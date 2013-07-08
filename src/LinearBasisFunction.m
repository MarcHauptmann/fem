classdef LinearBasisFunction
    properties
        left;
        center;
        right;
    end
    
    methods
        function self = LinearBasisFunction(left, center, right)            
            self.left = left;
            self.center = center;
            self.right = right;
        end
        
        function y = evaluate(self, x)
            y = self.evalForAll(@self.getValue, x);
        end
        
        function y = evalForAll(self, func, x)
            n = length(x);
            y = zeros(1, n);
            
            parfor i=1:n
                y(i) = func(x(i));
            end
        end
        
        function y = getValue(self, x)
            if x > self.left && x <= self.center
                y = (x-self.left)/(self.center-self.left);
            elseif x > self.center && x < self.right
                y = 1-(x-self.center)/(self.right-self.center);
            else
                y = 0;
            end
        end
        
        function y = evaluateDerivate(self, x)
            y = self.evalForAll(@self.getDerivate, x);
        end
        
        function y = getDerivate(self, x)
            if x > self.left && x <= self.center
                y = 1/(self.center - self.left);
            elseif x > self.center && x < self.right
                y = -1/(self.right - self.center);
            else
                y = 0;
            end
        end
    end
end

