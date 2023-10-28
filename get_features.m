function [features] = get_features(sig)
% Returns the extracted features from the signal
    
    fs = 4000;

    % the signal was cropped at first but decided to uncrop it due to
    % dimension error in cal_score function
    % cropping the start and the end of the signal
    % noise = 2*fs;
    % cropped_sig = sig(noise:end-noise);
    cropped_sig = sig;


    n = length(sig);
    
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
    
    normenv = normalize(filt2data);
   
    % getting stft 
    stft_sig = stft(normenv,fs,'frequencyrange', 'onesided', 'Window',kaiser(256,5),'OverlapLength', 220);
    full_stft = 0;
    for i=1:size(stft_sig,1)
      for j=1:size(stft_sig,2)
        full_stft = full_stft+abs(stft_sig(i,j));
      end
    end
    full_stft = fix(full_stft);

    % features from the frequency domain
    L = length(normenv); 
    Fn = fs/2;
    FT_af = fft(normenv)/L;
    Fv = linspace(0, 1, fix(L/2)+1)*Fn;
    Iv = 1:numel(Fv);
    [PksL,~] = findpeaks(abs(FT_af(Iv))*2);
    % max frequency peak
    max_peak_f = max(PksL);    

    % normalized frequency
    [pxx,~] = pwelch(normenv,'power');
    high_pwelch = max(pow2db(pxx));

    % wavelet
    [c,l] = wavedec(normenv,4,'db2');
    approx = appcoef(c,l,'db2');
    [cd1,cd2,cd3,cd4] = detcoef(c,l,[1 2 3 4]);

    entropy = mean(wentropy(logdata,'log energy'));

    % features being used
    features = [
                mean(normenv)
                median(normenv)
                mad(normenv)
                skewness(normenv)
                kurtosis(normenv)
                iqr(normenv)
                full_stft
                max_peak_f
                high_pwelch
                sum(abs(cd1))
                sum(abs(cd2))
                sum(abs(cd3))
                sum(abs(cd4))
                entropy
                ];

end




