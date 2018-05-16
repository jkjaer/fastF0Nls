function [harmonicDfts, pitchGrids] = ...
        computeComplexHarmonicDfts(dataVector, fullPitchGrid, ...
        pitchOrder, fftShiftVector)
    nDft = round(1/diff(fullPitchGrid(1:2)));
    dftData = fft(dataVector, nDft);
    % compensate for a symmetric time index (i.e., a start index of
    % -(nData-1)/2)
    fftShiftVectorLength = length(fftShiftVector);
    shiftedDftData = (dftData(1:fftShiftVectorLength).*fftShiftVector).';
    for ii = 1:pitchOrder
        dftIndices = computeDftIndicesNHarmonic(nDft, ...
            fullPitchGrid([1,end]), ii);
        nPitches = length(dftIndices);
        if ii == 1
            % allocate memory
            pitchGrids = nan(pitchOrder, nPitches);
            harmonicDfts = nan(pitchOrder, nPitches);
        end
        pitchGrids(ii,1:nPitches) = fullPitchGrid(1:nPitches);
        harmonicDfts(ii, 1:nPitches) = shiftedDftData(dftIndices+1);
    end
end

function dftIndices = computeDftIndicesNHarmonic(nDft, pitchBounds, ...
        pitchOrder)
    % the rounding operation is only used for numerical reasons in this
    % function
    minPitchIdx = max(0, round(pitchBounds(1)*nDft));
    % the max pitch must be smaller than 1/(2*modelOrder) to avoid
    % problems with aliasing
    maxPitchIdx = min(nDft/(2*pitchOrder)-1, ...
        round(pitchBounds(2)*nDft));
    dftIndices = (minPitchIdx:maxPitchIdx)*pitchOrder;
end
