classdef GenericFunctionTest < TestCase
    properties
    end
    
    methods
        function self = GenericFunctionTest(name)
            self = self@TestCase(name);
        end
        
        function testEvaluate_identity(self)
            f = GenericFunction(@(x) x);
            
            assertAlmostEqual(f.evaluate(1), 1);
            assertAlmostEqual(f.evaluate([1,2]), [1,2]);
        end
    end
end

