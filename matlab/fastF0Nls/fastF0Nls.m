classdef fastF0Nls < handle
%
% Computes an estimate of the fundamental frequency from a real-valued
% signal segment.
%
% The signal segment is modelled as
%  data(n) = a_0 + sum_i [ a_i cos(i 2 \pi f_0 n) - b_i sin(i 2 \pi f_0 n)] + e(n)
% 
% where the sum is over L harmonic components. The scalar a_0 is the 
% DC-value. By default this is not included in the model.
%     
% Constructor input:
%
%    N: Signal length (integer, scalar)    
%    
%    L: maximum model number of harmonics (i.e., order) that is expected
%    (positive integer, scalar)
%
%    pitchBounds: Lower and upper bounds on the fundamental frequency in
%           cycles/sample. The lower bound should not be set lower than
%           1/N and the upper bound can at most be 0.5 (i.e., the nyquist 
%           frequency). The bounds are specified as a 1x2 vector, e.g., 
%           [80/samplingFreq 400/samplingFreq].
%
%    DC: Optional (default false).
%          If true, includes a_0 (i.e., a DC-term) in the model
%
%    F: Optional (default 5*N*L)
%           FFT length (integer, scalar). Grid resolution is approximately
%           samplingFreq/F Hz. Do not change this unless you know what 
%           you are doing. If you need a finer resolution of the pitch
%           estimate, please set the tolerance of the 'estimate' method 
%           (see below).
%    
%    epsilon: regularization (default 0.0)
%           Adds epsilon*N to the diagonal of the linear
%           system. Set this to, e.g., 1e-4 if your lowest candidate 
%           fundamental frequency is smaller than 1 cycle/segment.
%    
% returns a fastF0Nls handle.
%
% Methods:
%
%   fastF0Nls.estimate(x, tolerance): 
%          Computes estimates of fundamental frequency and model order.
%       
%       Input:
%    
%           x: signal (vector, length N)    
%           tolerance: (Optional, default is 1/F). Refinement 
%               tolerance in cycles/sample.
%
%       Output:
%          
%           f0: (optional) estimated fundamental frequency in
%               cycles/sample.
%        
%           l: (optional) estimated number of harmonics/model order
% 
%           alpha: (optional) estimated linear parameters organised as
%               [a_0 a_1 ... a_l b_1 ... b_l]
%    
%
%   fastF0Nls.computeCostFunctions(x): Computes all the cost functions for
%   every candidate model order and fundamental frequency.
%
%     Input: 
%    
%           x: signal (vector, length N)    
%
%     Output:
%
%           costFunctionMatrix: a matrix with the objective, size
%               (nPitches x L). NaN indicates that the objective
%               has not been evaluated due to the given pitch_bounds or 
%               because it is only necessary to evalute 
%               f \in [pitchBounds(1), min(pitchBounds(2), 0.5/l)] for the
%               model orders l = 1,...,L.
%
%           fullPitchGrid: (optional) values for the frequency grid 
%               (real, size nPitches x 1)
%
% Based on the papers:
%
%      Fast fundamental frequency estimation: Making a statistically 
%      efficient estimator computationally efficient. Nielsen, Jesper Kjær;
%      Jensen, Tobias Lindstrøm; Jensen, Jesper Rindom; Christensen, 
%      Mads Græsbøll; Jensen, Søren Holdt. In: Signal Processing, 135, 
%      2017, pp. 188-197.
%
%      and
%
%      Bayesian Model Comparison With the g-Prior. Nielsen, Jesper Kjær; 
%      Christensen, Mads Græsbøll; Cemgil, Ali Taylan; Jensen, Søren 
%      Holdt. In: IEEE Transactions on Signal Processing, 62 (1), 2014, 
%      pp. 225-238.
%
% Example (analyse speech):
%
% % load the mono speech signal
% [speechSignal, samplingFreq] = audioread('roy.wav');
% nData = length(speechSignal);
%
% % set up
% segmentTime = 0.025; % seconds
% segmentLength = round(segmentTime*samplingFreq); % samples
% nSegments = floor(nData/segmentLength);
% f0Bounds = [80, 400]/samplingFreq; % cycles/sample
% f0Res = 0.1/samplingFreq; % cycles/sample
% maxNoHarmonics = 15;
% f0Estimator = fastF0Nls(segmentLength, maxNoHarmonics, f0Bounds);
%
% % do the analysis
% idx = 1:segmentLength;
% f0Estimates = nan(1,nSegments); % cycles/sample
% for ii = 1:nSegments
%     speechSegment = speechSignal(idx);
%     f0Estimates(ii) = f0Estimator.estimate(speechSegment, f0Res);
%     idx = idx + segmentLength;
% end
%
% Authors:  
%      Jesper Kjær Nielsen
%      Tobias Lindstrøm Jensen
% 
% This file was obtained from https://github.com/jkjaer/fastF0Nls. Please
% visit this page for more information on, e.g., licensing.
%
    
    % Can only be set in the constructor
    properties (SetAccess=immutable) 
        N
        L
        F
        pitchBoundsOuter = [0.0 0.5] % The outer bounds that are acceptable
        
        % default values
        epsilon = 0.0
        dcIsIncluded = false
        epsilon_ref = 0.0
        
        % Precomputed quantities
        crossCorrelationVectors
        Gamma1
        Gamma2
        fftShiftVector
    end

    properties (SetAccess=private)
        pitchBounds  % The current active bounds 
        fullPitchGrid
        validFftIndices
        defaultRefinementTol
        refinementTol
        gPriorParam = 3
        logPitchPdfs % the pitch pdf (column) for every candidate model order
        logModelPmf % an (L+1)x1 vector
        
        % computed point estimates
        estPitch % cycles/sample
        estOrder
    end

    methods
        
        function obj = fastF0Nls(N, L, pitchBounds, varargin)

            % validate input 
            if length(varargin) >= 1
                if islogical(varargin{1})
                    obj.dcIsIncluded = varargin{1};
                else
                    error('Argument 4 is not of type logical (true/false)')
                end
            end
            if length(varargin) >= 2
                if isscalar(varargin{2})
                    obj.F = varargin{2};
                else
                    error('Argument 5 is not a scalar')
                end
            else
                obj.F = 5*N*L;
            end
            obj.defaultRefinementTol = 1/obj.F;
            
            if length(varargin) >= 3
                if isscalar(varargin{3})
                    obj.epsilon = varargin{3};
                else
                    error('Argument 6 is not a scalar')
                end
            end
            
            if length(varargin) > 4
                error('Too many input arguments')
            end
            
            if ~isscalar(N)
                error('Input argument N is not a scalar')
            else
                obj.N = N;
            end
            if ~isscalar(L)
                error('Input argument L is not a scalar')
            else
                obj.L = L;
            end
            
            % check the pitch bounds
            if ~ (isvector(pitchBounds) && length(pitchBounds) == 2) 
                error(['Input argument pitchBound must be a 2-' ...
                       'vector.']);
            elseif pitchBounds(1) < obj.pitchBoundsOuter(1) || ...
                pitchBounds(2) > obj.pitchBoundsOuter(2)
                
                error(['Input argument pitchBounds must be within the ' ...
                       'bounds specified in the constructor (at ' ...
                       'least [0.0 0.5]).']);
            elseif pitchBounds(2) <= pitchBounds(1)
                error(['The upper pitch bound must be bigger than the', ...
                    ' bigger than the lower pitch bound.']);
            elseif pitchBounds(1)*L >= obj.pitchBoundsOuter(2)
                error(['The lower pitch bound or the maximum model', ...
                    ' order is too big. Their product must be smaller', ...
                    ' than 1/2.']);
            else
                obj.pitchBounds = pitchBounds;
            end
            
            if pitchBounds(1) < 1/N
                warning(['The lower pitch bound is set lower than one '...
                    ' period/segment. Inaccurate results might be '...
                    'produced - especially if you do not set the'...
                    ' regularisation parameter.']);
            end

            % Init
            F = obj.F;
            minFftIndex = ceil(F*pitchBounds(1));
            maxFftIndex = floor(F*pitchBounds(2));
            obj.validFftIndices = (minFftIndex:maxFftIndex)';
            obj.fullPitchGrid = obj.validFftIndices/F;
            nPitches = length(obj.fullPitchGrid);
            
            % cross-correlation vectors
            obj.crossCorrelationVectors = ...
                [N*ones(1, nPitches)/2 + N*obj.epsilon;...
                sin(pi*(1:2*L)'*obj.fullPitchGrid'*N)./...
                (2*sin(pi*(1:2*L)'*obj.fullPitchGrid'))];

            obj.fftShiftVector = ...
                exp(1i*2*pi*(0:ceil(F/2)-1)'*(N-1)/(2*F));

            % Compute Gamma for the T+H and T-H systems
            [obj.Gamma1, obj.Gamma2] = computeGamma(L, F, pitchBounds, ...
                obj.crossCorrelationVectors, nPitches, ...
                obj.validFftIndices, ...
                obj.dcIsIncluded);
        end
        
        function varargout = computeCostFunctions(obj, x)
            % Compute and returns the cost functions for l=1,...,L

            % validate input 
            if ~isvector(x) && length(x) == obj.N
                error(['First argument x must be vector of ' ...
                               'length N=', num2str(obj.N)]);
            end
            x = reshape(x, obj.N, 1); % make sure its a column
                                      % vector

                
            varargout{1} = computeAllCostFunctions(x, obj.L, ...
            	obj.fullPitchGrid, obj.fftShiftVector, ...
                obj.crossCorrelationVectors, obj.Gamma1, obj.Gamma2, ...
                obj.dcIsIncluded);
            
            if nargout == 2
                varargout{2} = obj.fullPitchGrid;
            end
            
        end
        
        function varargout = estimate(obj, x, varargin)
            % Estimate fundamental frequency and order of the signal x
                
            % validate input 
            if ~isvector(x) && length(x) == obj.N
                error(['First argument x must be vector of ' ...
                               'length N=', num2str(obj.N)]);
            end
            x = reshape(x, obj.N, 1); % make sure its a column
                                      % vector
            
            if length(varargin) == 1
                if isscalar(varargin{1})
                    obj.refinementTol = varargin{1};
                else
                    error('Argument 2 is not a scalar')
                end
            else
                obj.refinementTol = obj.defaultRefinementTol;
            end
            % if the data consists of all zeros then return zero model 
            % order
            if x'*x < 1e-14
                estimatedPitch = nan;
                estimatedOrder = 0;
            else
                % Step 1: compute the priors on the pitch and the model 
                % order for the current frame
                pitchLogPrior = -log(diff(obj.pitchBounds));
                logModelPrior = log(1/(obj.L+1));
                
                % Step 2: compute the profile log-likelihood function 
                % over only the fundamental frequency (the linear
                % parameters, noise variance and g have been integrated
                % out)
                costs = obj.computeCostFunctions(x);
                cod = costs*(1/(x'*x));
                [~, pitchLogLikelihood] = ...
                    computePitchLogLikelihood(cod, obj.N, obj.gPriorParam);
                
                % Step 3: compute the posteriors on the pitch and the model 
                % order for the current frame
                scaledPitchLogPosterior = ...
                    pitchLogLikelihood + pitchLogPrior;
                logMarginalLikelihood = computeLogMarginalLikelihood(...
                    scaledPitchLogPosterior, obj.fullPitchGrid);
                % normalise the pitch pdf
                obj.logPitchPdfs = scaledPitchLogPosterior-...
                    logMarginalLikelihood'*...
                    ones(1,size(scaledPitchLogPosterior,2));
                % compute the posterior model probabilities
                logMarginalLikelihood = [0, logMarginalLikelihood]; % add the null model
                postModelPmf = computePostModelPmf(...
                    logMarginalLikelihood, logModelPrior);
                obj.logModelPmf = log(postModelPmf);

                % step 4: compute point estimates of the model parameters
                [~, estimatedOrderIdx] = max(postModelPmf);
                estimatedOrder = estimatedOrderIdx-1;

                % Refine the pitch estimate if estimated model order > 0
                % and the tolerance is set to a smaller tolerance than what
                % is provided by the grid. Otherwise return the coarse
                % estimate.
                if estimatedOrder > 0
                    [~, pitchIndex] = ...
                        max(scaledPitchLogPosterior(estimatedOrder, :));
                    coarsePitchEstimate = obj.fullPitchGrid(pitchIndex(1));
                    if obj.refinementTol < obj.defaultRefinementTol
                        pitchLimits = ...
                            coarsePitchEstimate+[-1/obj.F, 1/obj.F];
                        costFunction = @(f0) -objFunction(f0, x, ...
                            estimatedOrder, obj.dcIsIncluded, ...
                            obj.epsilon_ref);
                        [f0l, f0u] = goldenSectionSearch(costFunction, ...
                            pitchLimits(1), pitchLimits(2), ...
                            obj.refinementTol);
                        estimatedPitch = (f0l+f0u)/2;
                    else
                        estimatedPitch = coarsePitchEstimate;
                    end
                else
                    estimatedPitch = nan;
                end
            end
            % store the computed estimates
            obj.estPitch = estimatedPitch;
            obj.estOrder = estimatedOrder;
            
            if nargout >= 1
                varargout{1} = estimatedPitch;
            end
            if nargout >= 2
                varargout{2} = estimatedOrder;
            end
            if nargout >= 3
                [~, alpha] = objFunction(estimatedPitch, x, ...
                    estimatedOrder, obj.dcIsIncluded);
                varargout{3} = alpha;
            end
            
        end
    end
end
