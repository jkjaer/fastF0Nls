function [marginalLikelihood, logMarginalLikelihood] = ...
        pitchMarginalLikelihoodLaplaceG(cod, nData, delta)
    [maxModelOrder, nFreqs] = size(cod);
    logMarginalLikelihood = nan(maxModelOrder, nFreqs);
    for iOrder = 1:maxModelOrder
        [gHat, tauVar] = computeLaplaceParameters(cod(iOrder,:), 1,  ...
            (nData-2*iOrder-delta)/2, nData/2);
        logMarginalLikelihood(iOrder,:) = ...
            log(gHat*(delta-2)/2)+(nData-2*iOrder-delta)/2*log(1+gHat)-...
            nData/2*log(1+gHat.*(1-cod(iOrder,:)))+...
            1/2*log(2*pi*tauVar);
    end
    marginalLikelihood = exp(logMarginalLikelihood);
end

function [gHat, tauVar] = computeLaplaceParameters(cod, v, w, u)
    a = (1-cod)*(v+w-u);
    b = (u-v)*cod+2*v+w-u;
    gHat = (b+sqrt(b.^2-4*a*v))./(-2*a);
    tauVar = 1./(gHat.*(1-cod)*u./(1+gHat.*(1-cod)).^2-...
        gHat.*w./(1+gHat).^2);
end

