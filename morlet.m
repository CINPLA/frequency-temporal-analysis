function [ Y, f_analysis ] = morlet(x, Fs, min_freq, max_freq, nFreq, K, varargin)
%MORLET(x, Fs, min_freq, max_freq, nFreq, K, varargin) morlet analysis
%
%
%   Performs morlet wavelet based bamd-pass filtering on x
%   x is a signal measured from t=0 to t=1. 
%   minfreq is the minimum frequency in the analysis
%   maxfreq is the maximum frequency in the analysis
%   nFreq is the number of frequencies in the analysis
%   K is the quality factor
%   the analysis does the convolution in the frequency domain and is of 
%   order n*log(n)*nFreq 
%   the filter frequencies are spread logarithmically, and the bin size
%   increases linearly as a function of analysis frequency.
p = inputParser;
def_logfreq = 0;
p.addRequired('x', @(x) isnumeric(x));
p.addRequired('Fs', @(x) isnumeric(x));
p.addRequired('minFreq', @(x) isnumeric(x));
p.addRequired('maxFreq', @(x) isnumeric(x));
p.addRequired('nFreq', @(x) isnumeric(x));
p.addRequired('K', @(x) isnumeric(x));
p.addParameter('logfreq', def_logfreq, @(x) isnumeric(x));
p.parse(x, Fs, min_freq, max_freq, nFreq, K, varargin{:});
logfreq = p.Results.logfreq;

n = length(x);

F_x = fft(x);

if logfreq
    f_analysis = logspace(log10(min_freq), log10(max_freq), nFreq); %analysis freqs
else
    f_analysis = linspace(min_freq, max_freq, nFreq); %analysis freqs
end

omega = linspace(0,Fs,length(x));


Y = zeros(n, nFreq); 
for i=1:nFreq
    f_a = f_analysis(i);
    c = sqrt(1+exp(-f_a*f_a) - 2*exp(-0.75*f_a*f_a));
%     c1 = sqrt(sqrt(2*pi)*f_a/K);
%     c2 = 1/c1;
    %k = exp(-0.5*sigma*sigma);
    F_psi = c*exp(-0.5*(f_a - omega).^2*K/f_a);

    y = F_x.*F_psi.'; 
    %size(y)
    %max(y)
    Y(:,i) = abs(ifft(y)); 
end




end

