%% Single-tone ENF estimation. Values of IFs bounded by [fc-bound, fc+bound]
%% Use this function only if input signal is properly bandpass filtered.
%   input:
%           signal:         time domain signal row vector
%           fs:             sampling frequency
%           window_dur:     window duration in second
%           step_size_dur:  step size duration in second
%           fc:             nominal single tone frequency (50 Hz, 100 Hz, etc.)
%           bound:          tolerable deviation of IF from fc in Hz
%           FFT_res_factor: number of FFT points = FFT_res_factor*fs
%   output:
%           IF0:            estimated ENF with interpolation
%           IF1:            estimated ENF wihtout interpolation
%           STFT_TFD:       estimated time-frequency distribution (TFD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IF1,IF0,STFT_TFD] = func_STFT_single_tone(signal,fs,window_dur,step_size_dur,fc,bound,FFT_res_factor)
window_length   = window_dur*fs;
window_func     = rectwin(window_length)';
step_size       = step_size_dur*fs;
NFFT            = FFT_res_factor*fs;
% signal          = [signal, zeros(1, window_length-1)]; % zero-padding
window_pos      = 1:step_size:(length(signal)-window_length+1);
STFT_TFD        = zeros (NFFT,length(window_pos));
IF0             = zeros(1,length(window_pos)); % output IF with interpolation
IF1             = zeros(1,length(window_pos)); % output IF without interpolation
absHalfTempFFT0 = zeros(1,NFFT/2);
N_in1           = round(NFFT/fs);
N_in2           = round(0.25*NFFT/fs);
for i = 1:length(window_pos)
    temp            = fft(signal(window_pos(i):window_pos(i)+window_length-1).*window_func,NFFT);
    STFT_TFD(:,i)   = temp;
    HalfTempFFT     = temp(1:end/2);
    absHalfTempFFT  = abs(HalfTempFFT);
    absHalfTempFFT0(fc*N_in1-N_in2:fc*N_in1+N_in2) = absHalfTempFFT(fc*N_in1-N_in2:fc*N_in1+N_in2);
    ValueMax        = max(absHalfTempFFT0);
    PeakLoc         = find(absHalfTempFFT0==ValueMax(1));
    PeakLoc         = PeakLoc(1);
    IF1(i)          = (PeakLoc-1)*fs/NFFT;
    ValueLeft       = HalfTempFFT(PeakLoc - 1);
    ValueCenter     = HalfTempFFT(PeakLoc);
    ValueRight      = HalfTempFFT(PeakLoc + 1);
    CorrectionCoef  = -real((ValueRight-ValueLeft)/(2*ValueCenter-ValueRight-ValueLeft));
    IF0(i)          = (PeakLoc+CorrectionCoef-1)*fs/NFFT; % interpolated estimation
end
IF1(IF1<fc-bound)=fc-bound;IF1(IF1>fc+bound)=fc+bound;
IF0(IF0<fc-bound)=fc-bound;IF0(IF0>fc+bound)=fc+bound;
end

