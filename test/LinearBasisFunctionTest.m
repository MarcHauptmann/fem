classdef LinearBasisFunctionTest < TestCase
    properties
        f;
        f2;
    end
    
    methods
        function self = LinearBasisFunctionTest(name)
            self = self@TestCase(name);
        end
        
        function setUp(self)
            self.f = LinearBasisFunction(0, 1, 2);
            self.f2 = LinearBasisFunction(0, 2, 4);
        end
        
        %% Tests für Funktionsauswertung
        
        function testBorder_intervalLengthIsOne_valuesAreZero(self)
            assertElementsAlmostEqual(self.f.evaluate(0), 0);
            assertElementsAlmostEqual(self.f.evaluate(2), 0);
        end
        
        function testCenter_intervalLengthIsOne_valueIsOne(self)
            assertElementsAlmostEqual(self.f.evaluate(1), 1);
        end
        
        function testLeftInterval_intervalLengthIsOne_valuesIncrease(self)
            assertElementsAlmostEqual(self.f.evaluate(0.3), 0.3);
            assertElementsAlmostEqual(self.f.evaluate(0.7), 0.7);
        end
        
        function testLeftInterval_intervalLengthIsTwo_valuesIncrease(self)
            assertElementsAlmostEqual(self.f2.evaluate(0.3), 0.15);
            assertElementsAlmostEqual(self.f2.evaluate(1.6), 0.8);
        end
        
        function testRightInterval_intervalLengthIsOne_valuesDecreases(self)
            assertElementsAlmostEqual(self.f.evaluate(1.3), 0.7);
            assertElementsAlmostEqual(self.f.evaluate(1.7), 0.3);
        end
        
        function testRightInterval_intervalLengthIsTwo_valuesDecreases(self)
            assertElementsAlmostEqual(self.f2.evaluate(2.4), 0.8);
            assertElementsAlmostEqual(self.f2.evaluate(3.6), 0.2);
        end
        
        function testVectorEvaluation_forValues(self)
            x = [-1, 0, 0.3, 0.7, 1, 1.3, 1.7, 2, 3];
            expected = [0, 0, 0.3, 0.7, 1, 0.7, 0.3, 0, 0];
            
            y = self.f.evaluate(x);
            
            assertElementsAlmostEqual(y, expected);
        end
        
        %% Tests für Auswertung der Ableitung
        
        function testLeftBorder_intervalLengthIsOne_derivatesAreZero(self)
            assertElementsAlmostEqual(self.f.evaluateDerivate(0), 0);
            assertElementsAlmostEqual(self.f.evaluateDerivate(2), 0);
        end
        
        function testLeftInterval_intervalLengthIsOne_derivateIsOne(self)
            assertElementsAlmostEqual(self.f.evaluateDerivate(0.3), 1);
            assertElementsAlmostEqual(self.f.evaluateDerivate(0.7), 1);
        end
        
        function testLeftInterval_intervalLengthIsTwo_derivateIsAHalf(self)
            assertElementsAlmostEqual(self.f2.evaluateDerivate(0.3), .5);
            assertElementsAlmostEqual(self.f2.evaluateDerivate(0.7), .5);
        end
        
        function testRightInterval_intervalLengthIsOne_derivateIsMinusOne(self)
            assertElementsAlmostEqual(self.f.evaluateDerivate(1.3), -1);
            assertElementsAlmostEqual(self.f.evaluateDerivate(1.7), -1);
        end
        
        function testRightInterval_intervalLengthIsTwo_derivateIsMinusAHalf(self)
            assertElementsAlmostEqual(self.f2.evaluateDerivate(2.3), -.5);
            assertElementsAlmostEqual(self.f2.evaluateDerivate(2.7), -.5);
        end
        
        function testVectorEvaluation_forDerivates(self)
            x = [-1, 1, 3, 5];
            expected = [0, 0.5, -0.5, 0];
            
            y = self.f2.evaluateDerivate(x);
            
            assertElementsAlmostEqual(y, expected);
        end
    end
end

