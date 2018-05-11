%% singlePitchMlCostFunctionNaiveImplementationMod
% A naive implementation of the single pitch ML cost-function.
%
%% Syntax:
%# [costFunctionMatrix,fullPitchGrid] = ...
%#      singlePitchMlCostFunctionNaiveImplementation(data, maxModelOrder, ...
%#      resolution, pitchBounds)
%
%% Description:
% A naive implementation of the single pitch ML cost-function.
% The cost function is computed on a uniform (Fourier) grid for all model
% orders up to the model orders specified by maxModelOrder. DC is also
% included in the model.
%
% * data: A data vector of real- or complex-valued data.
% * maxModelOrder: The maximum number of harmonic components that the cost
%   function is evaluated for.
% * resolution: The resolution of the uniform pitch grid that the cost functions
%   are evaluated for.
% * pitchBounds: A 2D vector specifiyng the lower and upper bound for the pitch
%   in cycles/sample.
% * costFunctionMatrix: A matrix of cost functions. Each row corresponds to the 
%   the cost function for a candidate model order.
% * fullPitchGrid: The pitches corresponding to the values in the rows of the 
%   cost function matrix.
%
%% Examples:
%# data = recordedData;
%# nData = length(data);
%# maxModelOrder = 10;
%# resolution = 1/(5*nData);
%# pitchBounds = [1/nData;0.5-1/nData];
%# [costFunctionMatrix,fullPitchGrid] = ...
%#     singlePitchMlCostFunctionNaiveImplementationMod(data, maxModelOrder, ...
%#     resolution, pitchBounds);
% 
function [costFunctionMatrix,fullPitchGrid, varargout] = ...
        singlePitchNLSCostsNaive(data,maxModelOrder,...
        resolution,pitchBounds,dcIsIncluded)
    
    nData = length(data);
    dataAreRealValued = isreal(data);
    time = (0:nData-1)'-(nData-1)/2;
    nFftGrid = ceil(1/resolution);
    minFftIndex = ceil(nFftGrid*pitchBounds(1));
    maxFftIndex = floor(nFftGrid*pitchBounds(2));
    validFftIndices = (minFftIndex:maxFftIndex)';
    fullPitchGrid = validFftIndices/nFftGrid;
    nPitches = length(fullPitchGrid);
    % initialisation
    costFunctionMatrix = nan(maxModelOrder,maxFftIndex-minFftIndex+1);
    if dcIsIncluded
        phiMatrix = costFunctionMatrix;
        dcMatrix = costFunctionMatrix;
    end
    
    for iOrder = 1:maxModelOrder
        if dataAreRealValued
            maxFftIndex = ...
                floor(min(nFftGrid*pitchBounds(2),nFftGrid/(2*iOrder)-1));
        else
            maxFftIndex = ...
                floor(min(nFftGrid*pitchBounds(2),nFftGrid/iOrder-1));
        end
        pitchGrid = (minFftIndex:maxFftIndex)'/nFftGrid;
        nPitches = length(pitchGrid);
        modelOrders = 1:iOrder;
        for jPitch = 1:nPitches
            pitch = pitchGrid(jPitch);
            expMatrix = exp(1i*2*pi*pitch*time*modelOrders);
            if dataAreRealValued
                if dcIsIncluded
                    sinusoidalMatrix = [ones(nData,1), real(expMatrix),...
                                        imag(expMatrix)];
                else
                    sinusoidalMatrix = [real(expMatrix),...
                                        imag(expMatrix)];
                end
            else
                if dcIsIncluded
                    sinusoidalMatrix = [ones(nData,1), expMatrix];
                else
                    sinusoidalMatrix = expMatrix;
                end
                
            end
            linearParameters = sinusoidalMatrix\data;
            if dcIsIncluded
                inverseFirstColumn = ...
                    (sinusoidalMatrix'*sinusoidalMatrix)\[1;zeros(2*iOrder,1)];
                dcMatrix(iOrder,jPitch) = linearParameters(1);
                phiMatrix(iOrder,jPitch) = inverseFirstColumn(1);
            end
            costFunctionMatrix(iOrder,jPitch) = ...
                real((data'*sinusoidalMatrix)*linearParameters);
        end
    end
    if dcIsIncluded
        varargout{1} = phiMatrix;
        varargout{2} = dcMatrix;
    else
        varargout = {};
    end
end
