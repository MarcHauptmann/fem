classdef createBasisTest < TestCase
    
    methods
        function self = createBasisTest(name)
            self = self@TestCase(name);
        end
        
        function testCreateBasis_linearBasis(self)
            % given
            elements = [0, 1, 2, 3];
            functions = {@(x) 1-x, @(x) x};
            x = linspace(0, 3, 7);
            
            B = [1 2
                 2 3
                 3 4];
            
            % when
            result = createBasis(functions, elements, x, B);
            
            % then
            expected = [  1   0   0   0
                        0.5 0.5   0   0
                          0   1   0   0
                          0 0.5 0.5   0
                          0   0   1   0
                          0   0 0.5 0.5
                          0   0   0   1];
                    
            assertEqual(result, expected);
        end
        
        function testCreateBasis_quadraticBasis(self)
            % given
            elements = [0, 1, 2, 3];
            f = {@(x) 1/2 * x.*(x-1)
                 @(x) 1-x.^2
                 @(x) 1/2 * x.*(x+1)};
            x = linspace(0, 3, 13);
            
            B = [1 2 3
                 3 4 5
                 5 6 7];
            
            % when
            result = createBasis(f, elements, x, B, -1, 1);
            
            % then
            i1 = [-1; -0.5; 0; 0.5; 1];
            i2 = [-0.5; 0; 0.5; 1];
            
            null2 = zeros(4, 1);
            null3 = zeros(5, 1);
            
            expected = [ f{1}(i1), f{2}(i1), f{3}(i1),    null3,    null3,    null3,    null3
                            null2,    null2, f{1}(i2), f{2}(i2), f{3}(i2),    null2,    null2
                            null2,    null2,    null2,    null2, f{1}(i2), f{2}(i2), f{3}(i2) ];
                           
            assertEqual(result, expected);
        end
    end
end

