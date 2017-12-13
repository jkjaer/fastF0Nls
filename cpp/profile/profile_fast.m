% A simple script to run the proposed fast algorithm
% The script generates data, start the profiler
% run the algorithm and show the profiler ourput for analysis
%
% Tobias L. Jensen, tlj@es.aau.dk
% June 2015
% Aalborg University, Denmark

clear all
addpath ../mlib
addpath ../sim/sim1

nData = 500;
maxModelOrder = 20;


modelOrder = 5;
modelOrders = (1:modelOrder)';
pitch = (1/modelOrder-0.5/nData)*rand(1)+0.5/nData;
pitchBounds = [1/nData,0.5-1/nData];
time0 = 0;
time = time0+(0:nData-1)';
Z = exp(1i*2*pi*pitch*time*modelOrders');
sinusoidalMatrix = [real(Z),-imag(Z)];
amplitudes = randn(2*modelOrder,1);
data = sinusoidalMatrix*amplitudes;
F = (nData*maxModelOrder*5);
resolution = 1/F;


profile on
w0 = methodWrapper(@singlePitchMlCostFunctionPrime,...
              data, maxModelOrder, resolution, pitchBounds, true);

%[actCostFunctionMatrix, fullPitchGrid] = ...
%    singlePitchMlCostFunctionPrime(data,...
%                              maxModelOrder,resolution,pitchBounds);

profile off
profile report