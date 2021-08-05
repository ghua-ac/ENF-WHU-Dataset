%% ENF harmonic estimation algorithms with harmonic enhancement
%   This is the program comparing multi-tone harmonic model based ENF estimation algorithms,
%   including
%       1. single-tone model for comparison
%       2. implementation of search within sum of harmonic components by Bykhovsky and Cohen [1]
%       3. implementation of search within weighted sum of harmonic component by Hajj-Ahmad et al. [2]
%       4. proposed method of harmonic selection
%
%   [1] D. Bykhovsky and A. Cohen, "Electrical network frequency (ENF) maximum-likelihood
%       estimation via a multitone harmonic model," IEEE Trans. Inf. Forensics Security,
%       vol. 8, no. 5, pp. 744�C753, May 2013.
%
%   [2] A. Hajj-Ahmad, R. Garg, and M. Wu, "Spectrum combining for ENF signal estimation,"
%       IEEE Signal Process. Lett., vol. 20, no. 9, pp. 885�C888, Sep. 2013.
%   For fair comparison no interpolation is used in frequency estimators.
%
%   ENF enhancement is achieved via robust filtering algorithm (RFA)
%   proposed by Hua and Zhang [3] for single-tone enhancement. Multi-tone extension of [3]
%   is provided in this program.
%
%   [3] G. Hua and H. Zhang, "ENF signal enhancement in audio recordings,"
%       IEEE Trans. Inf. Forensics Security, vol. 15, pp. 1868-1878, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
tic;
%% signal preparation
FS                       = 800; % constant sampling frequency
HARMONIC_INDEX           = [2,3,4,5,6,7]; % constant value for ENF harmonic processing
fc                       = 50*HARMONIC_INDEX; % nominal frequencies at each harmonic
bound                    = 0.1*HARMONIC_INDEX; % tolerable IF deviations at each harmonic
fundamental_f            = 50; % define ENF fundamental frequency
duration                 = 300; % second
N                        = FS*duration; % sample length
%     corrupted_index          = randperm(5,3)+2;
%     corrupted_index          = sort(corrupted_index); % indices of corrupted harmonics
corrupted_index          = [3, 6, 7];
% This is gonna take some time because each sample is generated according
% to all previous samples
[raw_wave_clean, IF_2nd] = func_ENF_synthesis_corrupted_harmonic(fundamental_f,HARMONIC_INDEX,corrupted_index,duration, FS);
SNR                      = -20; % dB
w                        = randn(1, N);
w                        = w/norm(w)/(10^(SNR/20));
raw_wave                 = raw_wave_clean + w;
disp(['Signal generated.']);
disp(['Corrupted harmonic indices:',num2str(corrupted_index),'.']);
filter_length            = 256;
[BPF_coeffs, coeffs_2nd] = func_BPF(filter_length);
input                    = filtfilt(BPF_coeffs,1,raw_wave);

%% harmonic enhancement
N_ite               = 2; % enhancement iterations
h_rfa               = 3000; % window length of RFA
initial_guess       = fc(1)*ones(1,N); % initial IF of RFA is fixed to 100 Hz
TS                  = 1; % constant time-step for RFA, 1 second
window_dur_rfa      = 8; % enhancement window size 8 seconds
FFT_res_rfa         = 200; % FFT resolution for RFA 1/FFT_res_rfa Hz
refined_guess       = initial_guess;
for k = 1:N_ite
    [input_denoised_single,~,refined_guess] = func_RFA(input,h_rfa,FS,TS,...
        refined_guess,fc(1),bound(1),window_dur_rfa,FFT_res_rfa);
end
initial_guesses     =  kron(HARMONIC_INDEX(2:end)'/2,refined_guess);
refined_guesses     = initial_guesses;
for k=1:N_ite
    [input_denoised_multi,~,refined_guesses] = func_RFA_multi(input,h_rfa,FS,TS,...
        refined_guesses,fc(2:end),bound(2:end),window_dur_rfa,FFT_res_rfa);
end
input_single_E      = input_denoised_single;
input_multi_E       = sum([input_denoised_single;input_denoised_multi],1);
input_original      = input;
disp(['Enhancement done.']);
%% ENF estimators
% set up parameters for frame-based processing
window_dur             = 16; % duration of overlapping frame in second
step_size_dur          = 1; % frame step-size usually 1 second
FFT_res_factor         = 2000; % FFT resolution = 1/FFT_res_factor Hz

f_ref                  = func_STFT_single_tone(raw_wave_clean,FS,window_dur,step_size_dur,fc(1),bound(1),FFT_res_factor);
correlation_thre       = func_empirical_coeff_thre(length(f_ref),10000)*4; % empirical correlation threshold
correlation_thre       = min(correlation_thre, 0.8);

f_single               = func_STFT_single_tone(input_original,FS,window_dur,step_size_dur,fc(1),bound(1),FFT_res_factor);
f_E_single             = func_STFT_single_tone(input_single_E,FS,window_dur,step_size_dur,fc(1),bound(1),FFT_res_factor);
f_MLE                  = func_STFT_multi_tone_search(input_original,FS,window_dur,step_size_dur,fc,bound,2*FFT_res_factor);
f_WMLE                 = func_STFT_multi_tone_search_weighted(input_original,FS,window_dur,step_size_dur,fc,bound,2*FFT_res_factor);

f_E_MLE                = func_STFT_multi_tone_search(input_multi_E,FS,window_dur,step_size_dur,fc,bound,2*FFT_res_factor);
f_E_WMLE               = func_STFT_multi_tone_search_weighted(input_multi_E,FS,window_dur,step_size_dur,fc,bound,2*FFT_res_factor);
[f_S_MLE, ~, s_index]  = func_STFT_multi_tone_MWC(input_original,FS,window_dur,step_size_dur,fc,bound,FFT_res_factor,correlation_thre,0);
f_S_WMLE               = func_STFT_multi_tone_MWC(input_original,FS,window_dur,step_size_dur,fc,bound,FFT_res_factor,correlation_thre,1);

[f_P_MLE, ~, p_index]  = func_STFT_multi_tone_MWC(input_multi_E,FS,window_dur,step_size_dur,fc,bound,FFT_res_factor,correlation_thre,0);
f_P_WMLE               = func_STFT_multi_tone_MWC(input_multi_E,FS,window_dur,step_size_dur,fc,bound,FFT_res_factor,correlation_thre,1);
toc;
MSE_single     = 1/length(f_ref)*norm(f_single-f_ref).^2
MSE_E_single   = 1/length(f_ref)*norm(f_E_single-f_ref).^2
MSE_MLE        = 1/length(f_ref)*norm(f_MLE-f_ref).^2
MSE_WMLE       = 1/length(f_ref)*norm(f_WMLE-f_ref).^2
MSE_E_MLE      = 1/length(f_ref)*norm(f_E_MLE-f_ref).^2
MSE_E_WMLE     = 1/length(f_ref)*norm(f_E_WMLE-f_ref).^2
MSE_S_MLE      = 1/length(f_ref)*norm(f_S_MLE-f_ref).^2
MSE_S_WMLE     = 1/length(f_ref)*norm(f_S_WMLE-f_ref).^2
MSE_P_MLE      = 1/length(f_ref)*norm(f_P_MLE-f_ref).^2
MSE_P_WMLE     = 1/length(f_ref)*norm(f_P_WMLE-f_ref).^2

%
% ENF_harmonics0  = func_STFT_multi_tone_component(input,FS,window_dur,step_size_dur,fc,bound,FFT_res_factor);
% ENF_harmonics1  = func_STFT_multi_tone_component(input_multi_E,FS,window_dur,step_size_dur,fc,bound,FFT_res_factor);
% NN              = length(f_ref);
% ll              = 49.9;
% hh              = 50.1;
% figure(2); plot(1:NN,ENF_harmonics0(1,:),1:NN,f_ref/2*2);axis([1 NN 2*ll 2*hh]);
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);grid on;
% figure(3);plot(1:NN,ENF_harmonics0(2,:),1:NN,f_ref/2*3);axis([1 NN 3*ll 3*hh]);grid on;
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
% figure(4);plot(1:NN,ENF_harmonics0(3,:),1:NN,f_ref/2*4);axis([1 NN 4*ll 4*hh]);grid on;
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
% figure(5);plot(1:NN,ENF_harmonics0(4,:),1:NN,f_ref/2*5);axis([1 NN 5*ll 5*hh]);grid on;
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
% figure(6);plot(1:NN,ENF_harmonics0(5,:),1:NN,f_ref/2*6);axis([1 NN 6*ll 6*hh]);grid on;
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
% figure(7);plot(1:NN,ENF_harmonics0(6,:),1:NN,f_ref/2*7);axis([1 NN 7*ll 7*hh]);grid on;
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
%
% figure(22);plot(1:NN,ENF_harmonics1(1,:),1:NN,f_ref/2*2);axis([1 NN 2*ll 2*hh]);grid on;
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
% figure(33);plot(1:NN,ENF_harmonics1(2,:),1:NN,f_ref/2*3);axis([1 NN 3*ll 3*hh]);grid on;
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
% figure(44);plot(1:NN,ENF_harmonics1(3,:),1:NN,f_ref/2*4);axis([1 NN 4*ll 4*hh]);grid on;
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
% figure(55);plot(1:NN,ENF_harmonics1(4,:),1:NN,f_ref/2*5);axis([1 NN 5*ll 5*hh]);grid on;
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
% figure(66);plot(1:NN,ENF_harmonics1(5,:),1:NN,f_ref/2*6);axis([1 NN 6*ll 6*hh]);grid on;
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
% figure(77);plot(1:NN,ENF_harmonics1(6,:),1:NN,f_ref/2*7);axis([1 NN 7*ll 7*hh]);grid on;
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
%
% figure(111); plot(1:NN,f_single,1:NN,f_ref/2*2);axis([1 NN 2*ll 2*hh]);
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);grid on;
% figure(112); plot(1:NN,f_E_single,1:NN,f_ref/2*2);axis([1 NN 2*ll 2*hh]);
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);grid on;
% figure(113); plot(1:NN,f_MLE,1:NN,f_ref/2*2);axis([1 NN 2*ll 2*hh]);
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);grid on;
% figure(114); plot(1:NN,f_WMLE,1:NN,f_ref/2*2);axis([1 NN 2*ll 2*hh]);
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);grid on;
% figure(115); plot(1:NN,f_E_MLE,1:NN,f_ref/2*2);axis([1 NN 2*ll 2*hh]);
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);grid on;
% figure(116); plot(1:NN,f_E_WMLE,1:NN,f_ref/2*2);axis([1 NN 2*ll 2*hh]);
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);grid on;
% figure(117); plot(1:NN,f_S_MLE,1:NN,f_ref/2*2);axis([1 NN 2*ll 2*hh]);
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);grid on;
% figure(118); plot(1:NN,f_S_WMLE,1:NN,f_ref/2*2);axis([1 NN 2*ll 2*hh]);
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);grid on;
% figure(119); plot(1:NN,f_P_MLE,1:NN,f_ref/2*2);axis([1 NN 2*ll 2*hh]);
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);grid on;
% figure(120); plot(1:NN,f_P_WMLE,1:NN,f_ref/2*2);axis([1 NN 2*ll 2*hh]);
% set(gca,'position',[0.01 0.01 0.98 0.98]);
% set(gca,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);grid on;

















