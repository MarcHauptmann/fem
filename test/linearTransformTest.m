classdef linearTransformTest < TestCase
    %LINEARTRANSFORMTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function self = linearTransformTest(name)
            self@TestCase(name);
        end
        
        function testLinearTransform_translation(self)
            % given
            x = [1, 1.5, 2];
            
            % when
            y = linearTransform(1, 2, 0, 1, x);
            
            % then
            expected = [0, 0.5, 1];
            
            assertEqual(y, expected);
        end
        
        function testLinearTransform_scaling(self)
            % given
            x = [1, 1.5, 2];
            
            % when
            y = linearTransform(1, 2, 1, 3, x);
            
            % then
            expected = [1, 2, 3];
            
            assertEqual(y, expected);
        end
        
        function testLinearTransform_fullTransformation(self)
            % given
            x = [1, 1.5, 2];
            
            % when
            y = linearTransform(1, 2, -1, 1, x);
            
            % then
            expected = [-1, 0, 1];
            
            assertEqual(y, expected);
        end
    end
    
end

