function [harmonicDfts, pitchGrids] = ...
        computeComplexHarmonicDfts(dataVector, fullPitchGrid, ...
        pitchOrder, fftShiftVector)
    nData = length(dataVector);
    dataAreRealValued = isreal(dataVector);
    nDft = round(1/diff(fullPitchGrid(1:2)));
    dftData = fft(dataVector, nDft);
    % compensate for a symmetric time index (i.e., a start index of
    % -(nData-1)/2)
    fftShiftVectorLength = length(fftShiftVector);
    shiftedDftData = (dftData(1:fftShiftVectorLength).*fftShiftVector).';
    for ii = 1:pitchOrder
        dftIndices = computeDftIndicesNHarmonic(nDft, ...
            fullPitchGrid([1,end]), ii, dataAreRealValued);
        nPitches = length(dftIndices);
        if ii == 1
            % allocate memory
            pitchGrids = nan(pitchOrder, nPitches);
            harmonicDfts = nan(pitchOrder, nPitches);
        end
        pitchGrids(ii,1:nPitches) = dftIndices/(ii*nDft);
        harmonicDfts(ii, 1:nPitches) = shiftedDftData(dftIndices+1);
    end
end

function dftIndices = computeDftIndicesNHarmonic(nDft, pitchBounds, ...
        pitchOrder, dataAreRealValued)
    minPitchIdx = max(0, ceil(pitchBounds(1)*nDft));
    % the max pitch must be smaller than 1/modelOrder for complex-valued
    % data and 1/(2*modelOrder) for real-valued data to avoid problems
    % with aliasing
    if dataAreRealValued
        maxPitchIdx = min(nDft/(2*pitchOrder)-1, ...
            floor(pitchBounds(2)*nDft));
    else
        maxPitchIdx = min(nDft/pitchOrder-1, ...
            floor(pitchBounds(2)*nDft));
    end    
    dftIndices = (minPitchIdx:maxPitchIdx)*pitchOrder;
end
