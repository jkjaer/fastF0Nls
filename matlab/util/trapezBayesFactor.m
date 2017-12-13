function bayesFactor = trapezBayesFactor(logMarginalLikelihood, ...
        fullPitchGrid)
    % compute the unnormalised Bayes' factors
    deltaFreq = diff(fullPitchGrid(1:2));
    maxModelOrder = size(logMarginalLikelihood,1);
    logBayesFactor = ones(1,maxModelOrder);
    for iOrder = 1:maxModelOrder
        % find the non-NaN indices
        indices = find(~isnan(logMarginalLikelihood(iOrder,:)));
        % find the maximum value
        maxLogMarginalLikelihood = max(logMarginalLikelihood(iOrder,indices));
        logBayesFactor(iOrder) = maxLogMarginalLikelihood+...
            log(deltaFreq*trapz(exp(...
            logMarginalLikelihood(iOrder,indices)-maxLogMarginalLikelihood)));
    end
    % add the null model
    logBayesFactor = [0,logBayesFactor];
    % normalise the Bayes' factors
    bayesFactor = exp(logBayesFactor-max(logBayesFactor));
    bayesFactor = bayesFactor/sum(bayesFactor);
end