clear
clc
close all

% add the estimator object to the MATLAB path
addpath ../fastF0Nls

% load the mono speech signal
[speechSignal, samplingFreq] = audioread('roy.wav');
nData = length(speechSignal);

% set up
segmentTime = 0.025; % seconds
segmentLength = round(segmentTime*samplingFreq); % samples
nSegments = floor(nData/segmentLength);
f0Bounds = [80, 400]/samplingFreq; % cycles/sample
maxNoHarmonics = 10;
f0Estimator = fastF0Nls(segmentLength, maxNoHarmonics, f0Bounds);

% do the analysis
idx = 1:segmentLength;
f0Estimates = nan(1,nSegments); % cycles/sample
for ii = 1:nSegments
    disp(['Processing segment ', num2str(ii), ' of ', ...
        num2str(nSegments)]);
    speechSegment = speechSignal(idx);
    f0Estimates(ii) = f0Estimator.estimate(speechSegment);
    idx = idx + segmentLength;
end
timeVector = (1:nSegments)*segmentTime-segmentTime/2;

%% compute the spectrogram of the signal
window = gausswin(segmentLength);
nOverlap = round(3*segmentLength/4);
nDft = 2048;
[stft, stftFreqVector, stftTimeVector] = ...
    spectrogram(speechSignal, window, nOverlap, nDft, samplingFreq);
powerSpectrum = abs(stft).^2;

%% plot the results
maxDynamicRange = 60; % dB
imagesc(stftTimeVector, stftFreqVector, ...
    10*log10(dynamicRangeLimiting(powerSpectrum, maxDynamicRange)));
set(gca,'YDir','normal')
hold on
plot(timeVector, f0Estimates*samplingFreq, 'r.')
hold off
colorbar
xlabel('time [s]')
ylabel('frequency [Hz]')
title('Why where you away a year, Roy?')