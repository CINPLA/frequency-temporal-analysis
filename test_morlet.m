close all; clear all;
[x,Fs] = wavread('fiolin.wav'); 
x = x(1:2^18,1);
n = length(x);
minFreq = 100; 
maxFreq = 2000; 
nFreq = 1000;
K = 50;


[Y, freqs] = morlet( x, Fs, minFreq, maxFreq, nFreq, K);

imagesc((1:n)/Fs, log10(freqs), Y')
set(gca, 'ydir', 'normal')
xlabel('time')
ylabel('log10(frequency)')


% cone of influence
R = K./freqs;
hold on
plot(R,log10(freqs),'LineWidth',2,'Color', 'k')
plot(n/Fs - R, log10(freqs),'LineWidth',2,'Color', 'k')
