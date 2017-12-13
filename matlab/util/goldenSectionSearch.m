%% goldenSectionSearch
% Single variable bounded minisation using a golden section search
%
%% Syntax:
% [lowerFinalBound, upperFinalBound] = goldenSectionSearch(objectiveFunction,...
%     lowerInitialBound, upperInitialBound, requiredInterval);
%
%% Description:
% Narrows down the minimum of a function from an initial interval to smaller
% interval using af golden section search.
% The implementation is based on Sec. 4.4 in 
% A. Antoniou and W.-S. Lu, Practical Optimization: Algorithms and Engineering
% Applications. Springer, Mar. 2007.
% * objectiveFunction: A function handle of the function to be minimised
% * lowerInitialBound: The lower initial bound of the minimiser
% * upperInitialBound: The upper initial bound of the minimiser
% * requiredInterval: The required final interval of uncertainty
% * lowerFinalBound: The final lower bound
% * upperFinalBound: The final upper bound
%% Examples:
% f = @(x) x.^2-2*x-3;
% lowerBound = -2;
% upperBound = 5;
% requiredUncertainty = 1e-3;
% [lowerBound, upperBound] = goldenSectionSearch(f,lowerBound,upperBound,...
%     requiredUncertainty);
% estimatedMinimiser = (upperBound+lowerBound)/2;
%
function [lowerFinalBound, upperFinalBound] = goldenSectionSearch(objectiveFunction,...
        lowerInitialBound, upperInitialBound, requiredInterval)
    % Check the input arguments
    if nargin<4
        error('goldenSectionSearch:argChk','Four input arguments are required.');
    end
    validateObjectiveFunction(objectiveFunction);
    validateBound(lowerInitialBound);
    validateBound(upperInitialBound);
    initialInterval = upperInitialBound-lowerInitialBound;
    if initialInterval<=0
        error('goldenSectionSearch:argChk',...
            'The lower bound must be smaller than the upper bound.');
    end
    validateInterval(requiredInterval);
    % Perform the golden section search
    [lowerFinalBound, upperFinalBound] = narrowBounds(objectiveFunction,...
        lowerInitialBound, upperInitialBound, requiredInterval);
end

% Return an error if the objective function is invalid
function validateObjectiveFunction(objectiveFunction)
    if ~isa(objectiveFunction, 'function_handle')
        error('goldenSectionSearch:argChk',...
                'The input argument must be a function handle');
    end
end

% Return an error if the bound is invalid
function validateBound(bound)
    if ~isscalar(bound) || ~isfinite(bound)
        error('goldenSectionSearch:argChk',...
            'The bounds must be real-valued scalars.');
    end
end

% Return an error if the interval is invalid
function validateInterval(interval)
    if ~isscalar(interval) || ~isfinite(interval) ||...
            interval<=0
        error('goldenSectionSearch:argChk',...
            'The required interval must be a positive real-valued scalar.');
    end
end

% Run the golden section search algorithm to narrow down the lower and upper 
% bound for the minimiser to within a tolerance of at most +/-requiredInterval/2. 
% The algorithm is based on algorithm 4.2 in A. Antoniou and W.-S. Lu, 
% Practical Optimization: Algorithms and Engineering Applications. 
% Springer, Mar. 2007.
function [lowerBound, upperBound] = narrowBounds(objectiveFunction,...
        lowerBound, upperBound, requiredInterval)
    goldenRatio = 0.5+sqrt(1.25);
    startInterval = upperBound-lowerBound;
    iInterval = startInterval/goldenRatio;
    variableLowerVal = upperBound-iInterval;
    funcLowerVal = objectiveFunction(variableLowerVal);
    variableUpperVal = lowerBound+iInterval;
    funcUpperVal = objectiveFunction(variableUpperVal);
    while iInterval > requiredInterval
        iInterval = iInterval/goldenRatio;
        if funcLowerVal > funcUpperVal
            % The minimum is in the interval [variableLowerVal;upperBound]
            lowerBound = variableLowerVal;
            variableLowerVal = variableUpperVal;
            variableUpperVal = lowerBound+iInterval;
            funcLowerVal = funcUpperVal;
            funcUpperVal = objectiveFunction(variableUpperVal);
        else
            % The minimum is in the interval [lowerBound;variableUpperVal]
            upperBound = variableUpperVal;
            variableUpperVal = variableLowerVal;
            variableLowerVal = upperBound-iInterval;
            funcUpperVal = funcLowerVal;
            funcLowerVal = objectiveFunction(variableLowerVal);
        end
        % If the final search tolerance is within the precision of
        % the computer, then stop
        if variableLowerVal>variableUpperVal
            break;
        end
    end
    % When iInterval is smaller than the required interval, we make one final
    % comparison so that (upperBound-lowerBound)<requiredInterval
    if funcLowerVal > funcUpperVal
        % The minimum is in the interval [variableLowerVal;upperBound]
        lowerBound = variableLowerVal;
    elseif funcLowerVal == funcUpperVal
        % The minimum is in the interval [variableLowerVal;variableUpperVal]
        lowerBound = variableLowerVal;
        upperBound = variableUpperVal;
    else
         % The minimum is in the interval [lowerBound;variableUpperVal]
        upperBound = variableUpperVal;
    end
end
