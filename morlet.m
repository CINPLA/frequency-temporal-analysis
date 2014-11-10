function [ Y, f_as ] = morlet( x, Fs,minFreq, maxFreq, nFreq, K)
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


n = length(x);

F_x = fft(x);

f_as = linspace(minFreq, maxFreq, nFreq); %analysis freqs

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

