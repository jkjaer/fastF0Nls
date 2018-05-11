function tests = computeComplexHarmonicDftsTest()
    tests = functiontests(localfunctions);
end

function testComputedDftsComplex(testCase)
    % setup
    rng(2);
    nData = 200;
    startIndex = -(nData-1)/2;
    dataVector = randn(nData,2)*[1;1i];
    maxPitchOrder = 5;
    pitchResolution = 1/(5*maxPitchOrder*nData);
    pitchBounds = [0,nData-1]/nData;
    fullPitchGrid = computePitchGrid(pitchResolution, pitchBounds, 1);
    nDft = round(1/diff(fullPitchGrid(1:2)));
    fftShiftVector = exp(-1i*2*pi*((0:nDft-1)/nDft)*startIndex);
    for ii = 1:maxPitchOrder
        pitchGrid = computePitchGrid(pitchResolution, pitchBounds, ii);
        nPitches = length(pitchGrid);
        if ii == 1
            expHarmonicDfts = nan(maxPitchOrder, nPitches);
            expPitchGrids = nan(maxPitchOrder, nPitches);
        end
        expPitchGrids(ii,1:nPitches) = pitchGrid;
        for jj = 1:nPitches
            sinusoid = ...
                exp(1i*2*pi*pitchGrid(jj)*ii*((0:nData-1)+startIndex)');
            expHarmonicDfts(ii,jj) = sinusoid'*dataVector;
        end
    end
    [actHarmonicDfts, actPitchGrids] = ...
        computeComplexHarmonicDfts(dataVector, fullPitchGrid, ...
        maxPitchOrder, fftShiftVector);
    testCase.assertEqual(actHarmonicDfts, expHarmonicDfts, ...
        'absTol', 1e-11);
    testCase.assertEqual(actPitchGrids, expPitchGrids, ...
        'absTol', 1e-11);
end

function testComputedDftsReal(testCase)
    % setup
    rng(2);
    nData = 200;
    startIndex = -(nData-1)/2;
    dataVector = randn(nData,1);
    maxPitchOrder = 5;
    pitchResolution = 1/(5*maxPitchOrder*nData);
    pitchBounds = [0,nData/2-1]/nData;
    fullPitchGrid = ...
        computePitchGrid(pitchResolution, pitchBounds, 1, true);
    nDft = round(1/diff(fullPitchGrid(1:2)));
    fftShiftVector = exp(-1i*2*pi*((0:nDft-1)/nDft)*startIndex);
    for ii = 1:maxPitchOrder
        pitchGrid = ...
            computePitchGrid(pitchResolution, pitchBounds, ii, true);
        nPitches = length(pitchGrid);
        if ii == 1
            expHarmonicDfts = nan(maxPitchOrder, nPitches);
            expPitchGrids = nan(maxPitchOrder, nPitches);
        end
        expPitchGrids(ii,1:nPitches) = pitchGrid;
        for jj = 1:nPitches
            sinusoid = ...
                exp(1i*2*pi*pitchGrid(jj)*ii*((0:nData-1)+startIndex)');
            expHarmonicDfts(ii,jj) = sinusoid'*dataVector;
        end
    end
    [actHarmonicDfts, actPitchGrids] = ...
        computeComplexHarmonicDfts(dataVector, fullPitchGrid, ...
        maxPitchOrder, fftShiftVector);
    testCase.assertEqual(actHarmonicDfts, expHarmonicDfts, ...
        'absTol', 1e-11);
    testCase.assertEqual(actPitchGrids, expPitchGrids, ...
        'absTol', 1e-11);
end

function pitchGrid = computePitchGrid(pitchResolution, pitchBounds, ...
        modelOrder, dataAreRealValued)
    if nargin < 4
        dataAreRealValued = false;
    end
    nPitches = ceil(1/pitchResolution);
    minPitchIdx = max(0, ceil(pitchBounds(1)*nPitches));
    % the max pitch must be smaller than 1/modelOrder to avoid problems
    % with aliasing
    if dataAreRealValued
        maxPitchIdx = min(nPitches/(2*modelOrder)-1, ...
            floor(pitchBounds(2)*nPitches));
    else
        maxPitchIdx = min(nPitches/modelOrder-1, ...
            floor(pitchBounds(2)*nPitches));
    end
    % the pitch grid must a subset of the full dft grid
    dftPitchGrid = (0:nPitches-1)'/nPitches;
    pitchGrid = dftPitchGrid((minPitchIdx(1):maxPitchIdx(1))+1);
end