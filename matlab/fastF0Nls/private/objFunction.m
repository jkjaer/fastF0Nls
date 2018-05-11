function varargout  = objFunction(f0, x, l, dcIsIncluded, varargin)
%
% Computes the objective of the function
%
% f(f0, l) = x^T Z_l(f0)(Z_l(f0)^T Z_l(f0))^-1 Z_l^T(f0) x
%
% Input:
%     omega:     fundamental frequency in radians 
%     x:         x length N
%     l:         model order
%     epsilon:   Regularization (default epsilon=0)
%
% Output:
%     objective: f(\omega, l)
%     alpha    : (optional) least-squares solution
%

N = length(x);

epsilon = 0;
if length(varargin) == 2
    epsilon = N*varargin{2};
end

n = (0:N-1)'-(N-1)/2;

E = exp(1i*2*pi*f0*n*(1:l));
if dcIsIncluded
    C = [ones(N, 1), real(E)]; %cos(omega*n'*(1:l));
else
    C = real(E);
end

S = imag(E); %sin(omega*n'*(1:l));
bc = C'*x;
bs = -S'*x;

ac = (C'*C+epsilon*eye(size(C, 2)))\bc; % just solve it using an O(n^3) algorithm for now
as = (S'*S+epsilon*eye(l))\bs; 

obj = bc'*ac + bs'*as;
varargout{1} = obj;

if nargout >= 2
    alpha = [ac; as];
    varargout{2} = alpha;
end

