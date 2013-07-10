classdef GenericFunction
    
    properties
        func;
    end
    
    methods
        function self = GenericFunction(func)
            self.func = func;
        end
        
        function y = evaluate(self, x)
            y = self.func(x);
        end
    end    
end

