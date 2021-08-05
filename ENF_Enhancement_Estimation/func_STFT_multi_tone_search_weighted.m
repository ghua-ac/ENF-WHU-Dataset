%% ENF estimation based on search within weighted sum of harmonic components
%
% This function is the implementation of the maximum-likelihood estimator
% proposed by Bykhovsky and Cohen [2]
%
%   [2] A. Hajj-Ahmad, R. Garg, and M. Wu, "Spectrum combining for ENF signal estimation,"
%       IEEE Signal Process. Lett., vol. 20, no. 9, pp. 885¨C888, Sep. 2013.
%
%   input:
%           signal:         time domain signal row vector
%           fs:             sampling frequency
%           window_dur:     window duration in second
%           step_size_dur:  step size duration in second
%           fc:             nominal multi-tone frequencies [100,150,200,...]
%           bound:          tolerable deviation of IF from fc in Hz
%           FFT_res_factor: number of FFT points = FFT_res_factor*fs
%   output:
%           IF:            estimated ENF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IF,weights] = func_STFT_multi_tone_search_weighted(signal,fs,window_dur,step_size_dur,fc,bound,FFT_res_factor)
window_length   = window_dur*fs;
window_func     = rectwin(window_length)';
step_size       = step_size_dur*fs;
NFFT            = FFT_res_factor*fs;
% signal          = [signal, zeros(1, window_length-1)]; % zero-padding
window_pos      = 1:step_size:(length(signal)-window_length+1);
IF              = zeros(1,length(window_pos)); % output IF without interpolation
%% set bandwidth to estimate local SNR according to [2]
band_low_signal_1st   = round(49.98*FFT_res_factor);
band_high_signal_1st  = round(50.02*FFT_res_factor); % signal band: 49.98~50.02 Hz
band_low_noise_1st    = round(49.9*FFT_res_factor);
band_high_noise_1st   = round(50.1*FFT_res_factor); % noise band: 49.9~50.1 Hz excluding signal band                   
weights               = zeros(1,length(fc));
%% set harmonic search region
search_region_1st = round((50-bound(1)/2)*FFT_res_factor):round((50+bound(1)/2)*FFT_res_factor); % fundamental IF search region
length_per_band   = length(search_region_1st);
search_region     = kron(search_region_1st,(fc/50))'; % harmonic IF search region index
%% search loop
for i = 1:length(window_pos)
    temp           = fft(signal(window_pos(i):window_pos(i)+window_length-1).*window_func,NFFT);
    HalfTempFFT    = temp(1:end/2);
    absHalfTempFFT = abs(HalfTempFFT).';
    % calculate weights for each harmonic component according to [2]
    for j = 1:length(fc)
       signal_component   = absHalfTempFFT(band_low_signal_1st*(fc(j)/50):band_high_signal_1st*(fc(j)/50));
       signal_plus_noise  = absHalfTempFFT(band_low_noise_1st*(fc(j)/50):band_high_noise_1st*(fc(j)/50));
       weights(j)         = norm(signal_component).^2/(norm(signal_plus_noise).^2-norm(signal_component).^2);
    end
    weights        = weights/norm(weights);
    fbin_candidate = absHalfTempFFT(search_region);
    fbin_candidate = diag(weights)*reshape(fbin_candidate,[length(fc),length_per_band]);
    weighted_fbin  = sum(fbin_candidate.^2,1); % harmonically weighted frequency bin energy 
    ValueMax       = max(weighted_fbin);
    PeakLoc        = search_region_1st((weighted_fbin==ValueMax(1)));
    IF(i)          = PeakLoc*fs/NFFT*2; % (one sample shift compensated, location no need to "- 1")
end
norm_fc            = 100;
norm_bound         = 100*bound(1)/fc(1);
IF(IF<norm_fc-norm_bound)=norm_fc-norm_bound;IF(IF>norm_fc+norm_bound)=norm_fc+norm_bound;
end


