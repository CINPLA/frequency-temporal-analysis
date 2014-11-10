function [ Y, f_as ] = morlet(x, Fs, minFreq, maxFreq, nFreq, K, varargin)
%MORLET Performs morlet wavelet based bamd-pass filtering on x
%   x is a signal measured from t=0 to t=1. 
%   minfreq is the minimum frequency in the analysis
%   maxfreq is the maximum frequency in the analysis
%   nFreq is the number of frequencies in the analysis
%   sigma is the quality factor
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
p.addParamValue('logfreq', def_logfreq, @(x) isnumeric(x));
p.parse(x, Fs, minFreq, maxFreq, nFreq, K, varargin{:});
logfreq = p.Results.logfreq;

n = length(x);

F_x = fft(x);

if logfreq
    f_as = logspace(log10(minFreq), log10(maxFreq), nFreq); %analysis freqs
else
    f_as = linspace(minFreq, maxFreq, nFreq); %analysis freqs
end

omega = linspace(0,Fs,length(x));


Y = zeros(n, nFreq); 
for i=1:nFreq
    f_a = f_as(i);
    c = sqrt(1+exp(-f_a*f_a) - 2*exp(-0.75*f_a*f_a));
    
    %k = exp(-0.5*sigma*sigma);
    F_psi = c.*pi^0.25.*(exp(-0.5*(f_a - omega).^2*K));

    y = F_x.*F_psi.'; 
    %size(y)
    %max(y)
    Y(:,i) = sqrt(abs(ifft(y))); 
end




end

