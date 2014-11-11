function [ YY, f_analysis] = slider_plot_varied_K(x, Fs, min_freq, max_freq, nFreq, K_list, varargin)
%SLIDER_PLOT_VARIED_K(x, Fs, minFreq, maxFreq, nFreq, K, varargin) 
%   SLIDER_PLOT_VARIED_K(x, Fs, minFreq, maxFreq, nFreq, K, varargin) 
%   Creates a slider plot of morlet wavelet analyses with
%   different K values
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
p.addParamValue('logfreq', def_logfreq, @(x) isnumeric(x));
p.parse(x, Fs, min_freq, max_freq, nFreq, K_list, varargin{:});
logfreq = p.Results.logfreq;

n = length(x);

if logfreq
    f_analysis = logspace(log10(min_freq), log10(max_freq), nFreq); %analysis freqs
else
    f_analysis = linspace(min_freq, max_freq, nFreq); %analysis freqs
end




YY = zeros(n,nFreq,length(K_list));
for k=1:length(K_list)
    K = K_list(k);
    [Y, f_a] = morlet(x, Fs, min_freq, max_freq, nFreq, K);
    YY(:,:,k) = Y;
end


g = imagesc(YY(:,:,1)'); % g is now a pointer to a plot
set(gca,'ydir', 'normal')
title(sprintf('K=%i', K_list(1)));
h = uicontrol('style','slider','units','pixel','position',[0 0 200 20]);
addlistener(h,'ActionEvent',@(hObject, event) slider_event_listener(hObject, event,YY,g, K_list));


end



function slider_event_listener(hObject,event,X,plot_pointer, K_list)

n = get(hObject,'Value');
index_list = 1:length(K_list);
n = max(round(n*length(index_list)),1);


set(plot_pointer,'CData',X(:,:,n)');
title(sprintf('K=%i', K_list(n)))
drawnow;


end