function tests = test()
%
% Test script for evaluating the single pitch NLS cost function
%
%  Can be executed as runtests('test')
%
% Author:
%    J. K. Nielsen, jkn@es.aau.dk
%    T. L. Jensen,  tlj@es.aau.dk
%    Aalborg University, 2016
%   
    tests = functiontests(localfunctions);
end    


function testProducedCostFunctionNoDC(testCase)

    % Seeding may be used to make debugging easier.
    rng(100);

    N = 100;
    ell = 2;
    L = 10;
    alpha = [0; randn(2*ell,1)];
    tolerance = 1e-8;    
    pitchBounds = [1/N,0.5-1/N];
    f0 = (pi+1.5)/N;
    
    F = (N*L*5);
    resolution = 1/F;
    x = generator(N, ell, f0, alpha);

    [expCosts, ~] = ...
        singlePitchNLSCostsNaive(x, ...
        L, resolution, pitchBounds, false);

    f = fastF0Nls(N, L, pitchBounds);
    
    actCosts = f.computeCostFunctions(x);

    testCase.assertEqual(actCosts, expCosts, ...
                         'absTol', 1e-10);
    
    [f0h, ellh, alphah] = f.estimate(x, tolerance);

    testCase.assertEqual(f0h, f0, 'absTol', 1e-8)
    testCase.assertEqual(ellh, ell, 'absTol', 1e-14)
    testCase.assertEqual(alphah, alpha(2:end), 'absTol', 1e-8)
    
    % Try to set a new tolerance and make sure everything works
    tolerance = 1e-9; 

    actCosts = f.computeCostFunctions(x);

    testCase.assertEqual(actCosts, expCosts, ...
                         'absTol', 1e-10);
    
    [f0h, ellh] = f.estimate(x, tolerance);

    testCase.assertEqual(f0h, f0, 'absTol', 1e-8)
    testCase.assertEqual(ellh, ell, 'absTol', 1e-14)

end

function testProducedCostFunctionDC(testCase)

    % Seeding may be used to make debugging easier.
    rng(100);

    N = 100;
    ell = 2;
    L = 10;
    alpha = [2.1; randn(2*ell,1)];
    tolerance = 1e-8;  
    pitchBounds = [1/N, 0.5-1/N];
    f0 = (pi+1.5)/N;
    
    F = (N*L*5);
    resolution = 1/F;
    x = generator(N, ell, f0, alpha);

    [expCosts, ~] = ...
        singlePitchNLSCostsNaive(x, ...
        L, resolution, pitchBounds, true);

    f = fastF0Nls(N, L, pitchBounds, true);
    
    actCosts = f.computeCostFunctions(x);

    testCase.assertEqual(actCosts, expCosts, ...
                         'absTol', 1e-10);
    
    [f0h, ellh, alphah] = f.estimate(x, tolerance);

    testCase.assertEqual(f0h, f0, 'absTol', 1e-8)
    testCase.assertEqual(ellh, ell, 'absTol', 1e-14)
    testCase.assertEqual(alphah, alpha, 'absTol', 1e-8)
    
end

function testProducedCostFunctionDCHighRes(testCase)

    % Seeding may be used to make debugging easier.
    rng(100);

    N = 100;
    ell = 2;
    L = 10;
    alpha = [2.1; randn(2*ell,1)];
    tolerance = 1e-8;
    pitchBounds = [1/N, 0.5-1/N];
    f0 = (pi+1.5)/N;
    
    F = 10*N*L;
    resolution = 1/F;
    x = generator(N, ell, f0, alpha);

    [expCosts, ~] = ...
        singlePitchNLSCostsNaive(x, ...
        L, resolution, pitchBounds, true);

    f = fastF0Nls(N, L, pitchBounds, true, F);
    
    actCosts = f.computeCostFunctions(x);

    testCase.assertEqual(actCosts, expCosts, ...
                         'absTol', 1e-10);
    
    [f0h, ellh] = f.estimate(x, tolerance);

    testCase.assertEqual(f0h, f0, 'absTol', 1e-8)
    testCase.assertEqual(ellh, ell, 'absTol', 1e-14)
    
end

function testProducedCostFunctionZeroModelOrder(testCase)

    % Seeding may be used to make debugging easier.
    rng(100);

    N = 100;
    expModelOrder = 0;
    L = 10;
    tolerance = 1e-8;
    pitchBounds = [1/N, 0.5-1/N];
    expPitch = nan;
    
    F = 10*N*L;
    x = randn(N,1);

    f = fastF0Nls(N, L, pitchBounds, true, F);

    [actPitch, actModelOrder] = f.estimate(x, tolerance);

    testCase.assertEqual(actPitch, expPitch)
    testCase.assertEqual(actModelOrder, expModelOrder)
    
end
