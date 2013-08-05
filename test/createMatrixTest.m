classdef createMatrixTest < TestCase
    methods
        function self = createMatrixTest(name)
            self = self@TestCase(name);
        end
        
        function testCreateMatrix_linear_4elements(self)
            % given
            k = [1 -1 
                -1  1];
            
            geom = [1 2
                    2 3
                    3 4
                    4 5];
            
            % when
            result = createMatrix(k, geom);
            
            % then
            expected = [ 1 -1  0  0  0
                        -1  2 -1  0  0
                         0 -1  2 -1  0
                         0  0 -1  2 -1 
                         0  0  0 -1  1];
            
            assertEqual(expected, result);
        end
        
        function testCreateMatrix_quadratic_2elements(self)
            % given
            k = [ 7 -8  1
                 -8 16 -8
                  1 -8  7];
            
            geom = [1 2 3
                    3 4 5];
            
            % when
            result = createMatrix(k, geom);
            
            % then
            expected = [ 7 -8  1  0  0
                        -8 16 -8  0  0
                         1 -8 14 -8  1
                         0  0 -8 16 -8
                         0  0  1 -8  7];
            
            assertEqual(expected, result);
        end
    end
end

