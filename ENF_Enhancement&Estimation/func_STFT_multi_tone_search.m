%% ENF estimation based on search within sum of harmonic components
%
% This function is the implementation of the maximum-likelihood estimator
% proposed by Bykhovsky and Cohen [1]
%
%   [1] D. Bykhovsky and A. Cohen, "Electrical network frequency (ENF) maximum-likelihood
%       estimation via a multitone harmonic model," IEEE Trans. Inf. Forensics Security,
%       vol. 8, no. 5, pp. 744¨C753, May 2013.
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
function IF = func_STFT_multi_tone_search(signal,fs,window_dur,step_size_dur,fc,bound,FFT_res_factor)
window_length   = window_dur*fs;
window_func     = rectwin(window_length)';
step_size       = step_size_dur*fs;
NFFT            = FFT_res_factor*fs;
% signal          = [signal, zeros(1, window_length-1)]; % zero-padding
window_pos      = 1:step_size:(length(signal)-window_length+1);
IF              = zeros(1,length(window_pos)); % output IF without interpolation
%% set harmonic search region
search_region_1st = round((50-bound(1)/2)*FFT_res_factor):round((50+bound(1)/2)*FFT_res_factor); % fundamental IF search region
length_per_band   = length(search_region_1st);
search_region     = kron(search_region_1st,(fc/50))'; % harmonic IF search region index
%% search loop
for i = 1:length(window_pos)
    temp           = fft(signal(window_pos(i):window_pos(i)+window_length-1).*window_func,NFFT);
    HalfTempFFT    = temp(1:end/2);
    absHalfTempFFT = abs(HalfTempFFT).';
    fbin_candidate = absHalfTempFFT(search_region);
    fbin_candidate = reshape(fbin_candidate,[length(fc),length_per_band]);
    weighted_fbin  = sum(fbin_candidate.^2,1); % harmonically weighted frequency bin energy 
    ValueMax       = max(weighted_fbin);
    PeakLoc        = search_region_1st((weighted_fbin==ValueMax(1)));
    IF(i)          = PeakLoc*fs/NFFT*2; % (one sample shift compensated, location no need to "- 1")
end
norm_fc            = 100;
norm_bound         = 100*bound(1)/fc(1);
IF(IF<norm_fc-norm_bound)=norm_fc-norm_bound;IF(IF>norm_fc+norm_bound)=norm_fc+norm_bound;
end


