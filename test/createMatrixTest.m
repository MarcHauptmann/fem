classdef createMatrixTest < TestCase
    methods
        function self = createMatrixTest(name)
            self = self@TestCase(name);
        end
        
        function testCreateMatrix_linearBasis_4elements(self)
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
        
        function testCreateMatrix_linearBasis_4elements_oneH(self)
            % given
            k = [1 -1 
                -1  1];
            
            geom = [1 2
                    2 3
                    3 4
                    4 5];
            
            h = 4;
                
            % when
            result = createMatrix(k, geom, h);
            
            % then
            expected = 4 * [ 1 -1  0  0  0
                            -1  2 -1  0  0
                             0 -1  2 -1  0
                             0  0 -1  2 -1 
                             0  0  0 -1  1];
            
            assertEqual(expected, result);
        end
                
        function testCreateMatrix_linearBasis_4elements_hVector(self)
            % given
            k = [1 -1 
                -1  1];
            
            geom = [1 2
                    2 3
                    3 4];
            
            h = [1, 2, 1];
                
            % when
            result = createMatrix(k, geom, h);
            
            % then
            expected = [ 1 -1  0  0
                        -1  3 -2  0
                         0 -2  3 -1 
                         0  0 -1  1];
            
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

