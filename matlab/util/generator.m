function x = generator(N, ell, f0, alpha)
%
% Generates a harmonic signal
%

modelOrders = (1:ell)';
n = (-(N-1)/2:-(N-1)/2+N-1)';
Z = exp(1i*2*pi*f0*n*modelOrders');
sinusoidalMatrix = [ones(N, 1), real(Z), -imag(Z)];
x = sinusoidalMatrix*alpha;
