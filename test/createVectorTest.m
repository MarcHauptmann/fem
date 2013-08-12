classdef createVectorTest < TestCase    
    methods
        function self = createVectorTest(name)
            self = self@TestCase(name);
        end
        
        function testCreateVector_linear_constantF(self)
            % given
            elements = linspace(0, 1, 3);
            functions = {@(x) 1-x, @(x) x};
            geom = [ 1 2
                     2 3 ];
                 
             f= @(x) 1;
                 
            % when
            result = createVector(functions, geom, elements, f);
            
            % then
            expected = 0.5/2 * [1; 2; 1];
            
            assertAlmostEqual(expected, result);
        end
        
        function testCreateVector_linear_fIsFunction(self)
            % given
            elements = linspace(0, 1, 3);
            functions = {@(x) 1-x, @(x) x};
            geom = [ 1 2
                     2 3 ];
                 
            f = @(x) x;
                 
            % when
            result = createVector(functions, geom, elements, f);
            
            % then
            expected = [0.0416667
                        0.25
                        0.208333];
            
            assertAlmostEqual(expected, result, 0.001);
        end

        function testCreateVector_quadratic_constantF(self)
            % given
            elements = linspace(0, 1, 4);
            functions = {@(x) 1/2 * x.*(x-1)
                         @(x) 1-x.^2
                         @(x) 1/2 * x.*(x+1)};
            geom = [ 1 2 3
                     3 4 5 ];

            f = @(x) 1;

            % when
            result = createVector(functions, geom, elements, f, -1, 1);

            % then
            expected = 0.333/6 * [1; 4; 2; 4; 1];

            assertAlmostEqual(expected, result, 0.001);
        end
    end
end

