classdef fastF0Nls < handle
%
% Computes the non linear least-squares cost function of the real
% harmonic model with additive noise
%
%  data(n) = a_0 + sum_i a_i cos(i 2 \pi f_0 n) - b_i sin(i 2 \pi f_0 n) + e(n)
% 
% where L is the maxModelOrder. I.e. the objective is
%
%  f_l(w) = sum_n (a_0 sum_i^l a_i cos(i 2 \pi f n) 
%                - b_i sin(i 2 \pi f n) - data(n) )^2
%
% From the above objective the cost J_l(w) = sum_n data(n)^2 -
% f_l(w) is returned in CostFunctionMatrix for l=1:MaxModelOrder.
%     
% Constructor input:
%
%    N: Signal length (integer, scalar)    
%    
%    L: maxModelOrder that is considered (integer, scalar).
%
%    pitchBounds: Additional bounds of the cost function if tighter
%           bounds are known for the fundamental frequency w_0. A two
%           vector e.g. [0.001 0.49]. 0.0 is DC, 0.5 is the Nyquist.
%
%    DC: Optional (default false). 
%          If true, includes dc, a_0 in the model
%
%    F: Optional (default 5*N*L)
%           FFT length. Grid resolution is approximately 2/F (integer, scalar)
%    
%    epsilon: regularization (default 0.0)
%           Adds epsilon*N to the diagonal of the linear
%           system. Useful for large L.
%    
% returns a fastF0Nls handle.
%
% Methods:
%
%   fastF0Nls.computeCostFunctions(x): Computes all the cost functions
%
%     Input: 
%    
%         x: signal (vector, length N)    
%
%     Output:
%
%         costFunctionMatrix: a matrix with the objective, size
%           (nPitches x L). NaN indicates that the objective
%           has not been evaluated due to the given pitch_bounds or because
%           it is only necessary to evalute 
%           f \in [pitchBounds(1), min(pitchBounds(2), 0.5/l)] at
%           model order l = 1,...,maxModelOrder.
%
%         fullPitchGrid: (optional) values for the frequency grid 
%                (real, size nPitches x 1)
%   
%
%   fastF0Nls.estimate(x, tolerance): 
%          Computes estimates of fundamental frequency and model order.
%       
%      Input:
%    
%         x:    signal (vector, length N)    
%         tolerance: (Optional, default 1e-8). Refinement tolerance
%                       unit periods/sample
%
%     Output:
%          
%        f0:    (optional) estimated fundamental frequency. Units periods/sample
%        
%        l :    (optional) estimated model order
% 
%        alpha: (optional) linear parameters
%    
% This object implements a fast approach of evaluating the
%  objective function on a uniform grid for all model order.
%
% Current version: 1.1.0 (2018-01-09)
%
% Based on the papers:
%
%      Fast Fundamental Frequency Estimation: Making a
%      Statistically Efficient Estimator Computationally Efficient
%      J. K. Nielsen, T. L. Jensen, R. R. Jensen,
%      M. G. Christensen, S. H. Jensen.
%      Submitted, 2016
%
%      and
%
%      Fast and Statistically Efficient Fundamental Frequency Estimation
%      J. K. Nielsen, T. L. Jensen, R. R. Jensen,
%      M. G. Christensen, S. H. Jensen.
%      ICASSP, 2016
%
% Authors:  
%      J. K. Nielsen, kjn@es.aau.dk
%      T. L. Jensen, tlj@es.aau.dk
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

    properties
        tolerance = 1e-8
    end
    
    properties (SetAccess=private)
        pitchBounds  % The current active bounds 
        fullPitchGrid
        validFftIndices
        nPitches
    end

    methods
        
        function setTolerance(obj, tolerance)
            obj.tolerance = tolerance;
        end
                
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

            if ~ (isvector(pitchBounds) && length(pitchBounds) == 2) 
                error(['Input argument pitchBound must be a 2-' ...
                       'vector'])
            elseif pitchBounds(1) < obj.pitchBoundsOuter(1) || ...
                pitchBounds(2) > obj.pitchBoundsOuter(2)
                
                error(['Input argument pitchBounds must be within the ' ...
                       'bounds specified in the constructor (at ' ...
                       'least [0.0 0.5])'])
            else
                obj.pitchBounds = pitchBounds;
            end

            % Init
            F = obj.F;
            minFftIndex = ceil(F*pitchBounds(1));
            maxFftIndex = floor(F*pitchBounds(2));
            obj.validFftIndices = (minFftIndex:maxFftIndex)';
            obj.fullPitchGrid = obj.validFftIndices/F;
            nPitches = length(obj.fullPitchGrid);
            obj.nPitches = nPitches;
            
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
                error(sprintf(['First argument x must be vector of ' ...
                               'length N=%d'], obj.N))
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
            % Estimate fundamental and order of the input signal x
                
            % validate input 
            if ~isvector(x) && length(x) == obj.N
                error(sprintf(['First argument x must be vector of ' ...
                               'length N=%d'], obj.N))
            end
            x = reshape(x, obj.N, 1); % make sure its a column
                                      % vector
            
            if length(varargin) == 1
                if isscalar(varargin{1})
                    tolerance = varargin{1};
                else
                    error('Argument 2 is not a scalar')
                end
            else
                tolerance = obj.tolerance; % Use default
            end
            
            % if the data consists of all zeros then return zero model order
            if x'*x < 1e-14
                estimatedPitch = nan;
                estimatedOrder = 0;
            else
                % Step 1: Compute cost function
                costs = obj.computeCostFunctions(x);
                    
                % Step 2: Estimate model order and fundamental frequency
                cod = costs*(1/(x'*x));
                delta = 3;
                [~, logMarginalLikelihood] = ...
                    pitchMarginalLikelihoodLaplaceG(cod, obj.N, delta);
                bayesFactor = trapezBayesFactor(logMarginalLikelihood, ...
                                                obj.fullPitchGrid);
                [~, estimatedOrderIdx] = max(bayesFactor);
                estimatedOrder = estimatedOrderIdx-1;

                % Step 3: Refine if estimated model order > 0
                if estimatedOrder > 0
                    [~, pitchIndex] = max(costs(estimatedOrder, :));
                    coarsePitchEstimate = obj.fullPitchGrid(pitchIndex(1));
                    pitchLimits = coarsePitchEstimate+[-1/obj.F, 1/obj.F];
                    costFunction = @(f0) -objFunction(f0, x, ...
                        estimatedOrder, obj.dcIsIncluded, obj.epsilon_ref);
                    [f0l, f0u] = ...
                        goldenSectionSearch(costFunction, ...
                        pitchLimits(1), pitchLimits(2), tolerance);
                    estimatedPitch = (f0l+f0u)/2;
                else
                    estimatedPitch = nan;
                end
            end
                
            if nargout >= 1
                varargout{1} = estimatedPitch;
            end
            if nargout >= 2
                varargout{2} = estimatedOrder;
            end
            if nargout >= 3
                [~, alpha] = objFunction(estimatedPitch, x, estimatedOrder, obj.dcIsIncluded);
                varargout{3} = alpha;
            end
            
        end
        

    end
end
