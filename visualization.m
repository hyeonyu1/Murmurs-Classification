close all; 
clear; clc;
%%
% this is just a script to visualize the signal individually 
sig = audioread("data/normal/85033_TV.wav");
fs = 4000;

n = length(sig);
t = linspace(0, n/fs, n);

cropped_sig = sig;
%butterworth filter design and application
d = designfilt('bandpassiir','FilterOrder',2, ...
    'HalfPowerFrequency1',25,'HalfPowerFrequency2',400, ...
    'SampleRate',fs);

filtdata = filtfilt(d, cropped_sig);

%log scale:
logdata = -1*((filtdata-mean(filtdata)).^2).*log((filtdata-mean(filtdata)).^2);

%low pass filter design and application
d2 = designfilt('lowpassiir','FilterOrder',2, ...
    'PassBandFrequency',75, 'SampleRate',fs);

filtlogdata = filtfilt(d2, logdata);
    filt2data = sqrt(exp(filtlogdata));

normenv = normalize(filt2data);figure(1)
[pks, lcs] = findpeaks(normenv, 'MinPeakDistance', fs*0.15, 'MinPeakHeight', max(normenv)*0.03);

% Getting binary mask for heart sound
s12 = (0.12*fs)/2;
locs = lcs;
binary_mask = ones(size(sig));
for i=1:size(locs,1)
    lower_end = max(1, locs(i)-s12);
    upper_end = min(locs(i)+s12, size(sig,1));
    binary_mask(lower_end:upper_end,1) = 0;
end
    
plot(t, sig, t(lcs), sig(lcs), '*'); 
hold on; 
plot(t, binary_mask, 'black', 'LineWidth', 1);
xlim([0 n/fs]);
ylim([-0.5 2]);
grid on; 
xlabel("time(sec)");
legend("Cropped Sginal", "Peaks", "S1-S2 vs Systole-Diastole");
title("Signal with peaks and binary mask");




