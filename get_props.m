function [properties] = get_props(sig, features, visualize)
% Returns the properties of the signal
% If visualization is 1, it will show a plot

    fs = 4000;    

    % the signal was cropped at first but decided to uncrop it due to
    % dimension error in cal_score function
    % cropping the start and the end of the signal
    % noise = 2*fs;
    % sig = sig(noise:end-noise);

    n = length(sig);
    
    t = linspace(0, n/fs, n);
    
    %butterworth filter design and application
    d = designfilt('bandpassiir','FilterOrder',2, ...
        'HalfPowerFrequency1',25,'HalfPowerFrequency2',400, ...
        'SampleRate',fs);
    
    filtdata = filtfilt(d, sig);
    
    %log scale:
    logdata = -1*((filtdata-mean(filtdata)).^2).*log((filtdata-mean(filtdata)).^2);
    
    %low pass filter design and application
    d2 = designfilt('lowpassiir','FilterOrder',2, ...
        'PassBandFrequency',75, 'SampleRate',fs);
    
    filtlogdata = filtfilt(d2, logdata);
    filt2data = sqrt(exp(filtlogdata));
    
    % normalization
    normenv = normalize(filt2data);
    
    % peak locations:
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
    
    pathology_classifier = load("pathology_classifier.mat");

    data = array2table(features);
    data.Properties.VariableNames(1) = {'Var2'};    % mean(normenv)
    data.Properties.VariableNames(2) = {'Var3'};    % median(normenv)
    data.Properties.VariableNames(3) = {'Var4'};    % mad(normenv)
    data.Properties.VariableNames(4) = {'Var5'};    % skewness(normenv)
    data.Properties.VariableNames(5) = {'Var6'};    % kurtosis(normenv)
    data.Properties.VariableNames(6) = {'Var7'};    % iqr(normenv)
    data.Properties.VariableNames(7) = {'Var8'};    % full_stft
    data.Properties.VariableNames(8) = {'Var9'};    % max_peak_f
    data.Properties.VariableNames(9) = {'Var10'};   % high_pwelch
    data.Properties.VariableNames(10) = {'Var11'};  % sum(abs(cd1))
    data.Properties.VariableNames(11) = {'Var12'};  % sum(abs(cd2))
    data.Properties.VariableNames(12) = {'Var13'};  % sum(abs(cd3))
    data.Properties.VariableNames(13) = {'Var14'};  % sum(abs(cd4))
    data.Properties.VariableNames(14) = {'Var15'};  % entropy

    pathology =  pathology_classifier.pathology_classifier.predictFcn(data);
    properties.S_loc = lcs;
    properties.HR = ((numel(pks)/2)*60)/t(1,end);
    properties.ib_seg = binary_mask;
    properties.pathology = pathology;
    properties.len = length(sig);

    % Showing plot (pcg_plot)
    if visualize
        figure(1)
        plot(t, sig, t(lcs), sig(lcs), '*'); 
        hold on; 
        plot(t, binary_mask, 'black', 'LineWidth', 1);
        xlim([0 n/fs]);
        ylim([-0.5 2]);
        grid on; 
        xlabel("time(sec)");
        legend("Cropped Sginal", "Peaks", "S1-S2 vs Systole-Diastole");
        title("Signal with peaks and binary mask");
    end

end

