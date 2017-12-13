function limitedNonnegativeData = ...
        dynamicRangeLimiting(nonnegativeData, maxRangeDb)
    logPowerSpectrum = 10*log10(nonnegativeData);
    limitedLogPowerSpectrum = max(logPowerSpectrum, ...
            max(max(logPowerSpectrum)-maxRangeDb));
    limitedNonnegativeData = 10.^(limitedLogPowerSpectrum/10);
end