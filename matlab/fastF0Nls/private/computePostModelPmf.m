function postModelPmf = computePostModelPmf(logMarginalLikelihood, ...
        logModelPrior)
    scaledLogPostModelPmf = logMarginalLikelihood + logModelPrior;
    scaledPostModelPmf = ...
        exp(scaledLogPostModelPmf-max(scaledLogPostModelPmf));
    postModelPmf = scaledPostModelPmf/sum(scaledPostModelPmf);
end
