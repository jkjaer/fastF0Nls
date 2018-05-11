clear
clc
close all

% add the estimator object to the MATLAB path
addpath ../fastF0Nls

% audio is from the EBU SQAM CD
[rawData, rawSamplingFreq] = audioread('09viola.flac');
channelNo = 1;
% resample the signal to get faster processing
samplingFreq = 8000;
musicSignal = resample(rawData(:,channelNo), samplingFreq, ...
	rawSamplingFreq);
nData = length(musicSignal);

% set up
segmentTime = 0.03; % seconds
segmentLength = round(segmentTime*samplingFreq); % samples
overlap = 75; % between segment as a percentage
nShift = round((1-overlap/100)*segmentLength); % samples
shiftTime = nShift/samplingFreq; % seconds
nSegments = ceil((nData-segmentLength+1)/nShift);
f0Bounds = [100, 1000]/samplingFreq; % cycles/sample
maxNoHarmonics = 20;
f0Estimator = fastF0Nls(segmentLength, maxNoHarmonics, f0Bounds);
f0Tolerance = 0.1/samplingFreq; % tolerance of the refinement in cycles/sample

% do the analysis
idx = 1:segmentLength;
f0Estimates = nan(nSegments,1); % cycles/sample
orderEstimates = nan(nSegments,1);
for ii = 1:nSegments
    disp(['Processing segment ', num2str(ii), ' of ', ...
        num2str(nSegments)]);
    musicSegment = musicSignal(idx);
    [f0Estimates(ii), orderEstimates(ii)] = ...
        f0Estimator.estimate(musicSegment, f0Tolerance);
    idx = idx + nShift;
end
timeVector = segmentTime/2+(1:nSegments)*shiftTime-shiftTime/2;

%% compute the spectrogram of the signal
window = gausswin(segmentLength);
nOverlap = round(3*segmentLength/4);
nDft = 2048;
[stft, stftFreqVector, stftTimeVector] = ...
    spectrogram(musicSignal, window, nOverlap, nDft, samplingFreq);
powerSpectrum = abs(stft).^2;

%% plot the results
maxDynamicRange = 60; % dB
figure(1)
imagesc(stftTimeVector, stftFreqVector, ...
    10*log10(dynamicRangeLimiting(powerSpectrum, maxDynamicRange)));
set(gca,'YDir','normal')
hold on
plot(timeVector, f0Estimates*samplingFreq, 'r.')
hold off
colorbar
xlabel('time [s]')
ylabel('frequency [Hz]')

figure(2)
plot(timeVector, orderEstimates, '.')
xlabel('time [s]')
ylabel('estimated number of harmonics [.]')
