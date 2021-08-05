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
clear; close; clc;
%% parameter setting and initialization
FS                       = 800; % constant sampling frequency
HARMONIC_INDEX           = [2,3,4,5,6,7]; % constant value for ENF harmonic processing
fc                       = 50*HARMONIC_INDEX; % nominal frequencies at each harmonic
bound                    = 0.1*HARMONIC_INDEX; % tolerable IF deviations at each harmonic
filter_length            = 256;
[BPF_coeffs, coeffs_2nd] = func_BPF(filter_length);
index                    = dir('C:\Users\GHUA\Desktop\ENF-WHU-Dataset/H1/*.wav');
index_ref                = dir('C:\Users\GHUA\Desktop\ENF-WHU-Dataset/H1_ref/*.wav');

f_ref                    = cell(1,length(index));

f_single                 = cell(1,length(index));
f_E_single               = cell(1,length(index));
f_MLE                    = cell(1,length(index));
f_WMLE                   = cell(1,length(index));
f_E_MLE                  = cell(1,length(index));
f_E_WMLE                 = cell(1,length(index));
f_S_MLE                  = cell(1,length(index));
f_S_WMLE                 = cell(1,length(index));
f_P_MLE                  = cell(1,length(index));
f_P_WMLE                  = cell(1,length(index));

MSE_single               = zeros(1,length(index));
MSE_E_single             = zeros(1,length(index));
MSE_MLE                  = zeros(1,length(index));
MSE_WMLE                 = zeros(1,length(index));
MSE_E_MLE                = zeros(1,length(index));
MSE_E_WMLE               = zeros(1,length(index));
MSE_S_MLE                = zeros(1,length(index));
MSE_S_WMLE               = zeros(1,length(index));
MSE_P_MLE                = zeros(1,length(index));
MSE_P_WMLE               = zeros(1,length(index));

selected_harmonic_index0 = cell(1,length(index));
selected_harmonic_index  = cell(1,length(index));
tic;
for i =1:length(index)
    i
    %% import audio recording and reference
    [audio, fs_audio] = audioread(strcat('C:\Users\GHUA\Desktop\ENF-WHU-Dataset/H1/',index(i).name));
    [ref, fs_ref]     = audioread(strcat('C:\Users\GHUA\Desktop\ENF-WHU-Dataset/H1_ref/',index_ref(i).name));
    ref               = ref';
    audio             = audio(:,1);
    audio             = audio';
    raw_wave          = resample(audio, FS, fs_audio);
    N                 = length(raw_wave);
    %% bandpass filtering
    input             = filtfilt(BPF_coeffs,1,raw_wave); % multi-tone input signal without enhancement
    %% harmonic enhancement
    N_ite             = 2; % enhancement iterations
    h_rfa             = 3000; % window length of RFA
    initial_guess     = fc(1)*ones(1,N); % initial IF of RFA is fixed to 100 Hz
    TS                = 1; % constant time-step for RFA, 1 second
    window_dur_rfa    = 8; % enhancement window size 8 seconds
    FFT_res_rfa       = 200; % FFT resolution for RFA 1/FFT_res_rfa Hz
    refined_guess     = initial_guess;
    % step 1: single-tone enhancement
    for k = 1:N_ite
        [input_denoised_single,~,refined_guess] = func_RFA(input,h_rfa,FS,TS,...
            refined_guess,fc(1),bound(1),window_dur_rfa,FFT_res_rfa);
    end
    % step 2: use refined_guess to construct initial guesses for harmonics
    initial_guesses   =  kron(HARMONIC_INDEX(2:end)'/2,refined_guess);
    refined_guesses   = initial_guesses;
    % step 3: harmonic enhancement
    for k=1:N_ite
        [input_denoised_multi,~,refined_guesses] = func_RFA_multi(input,h_rfa,FS,TS,...
            refined_guesses,fc(2:end),bound(2:end),window_dur_rfa,FFT_res_rfa);
    end
    input_E_single    = input_denoised_single;
    input_E_multi     = sum([input_denoised_single;input_denoised_multi],1);
%   audiowrite(strcat('C:\Users\GHUA\Desktop\ENF-Forensic_Audio_Recording_Dataset/H1_enhanced/',index(i).name),input_multi,fs_audio);
    %% ENF estimators 
    % set up parameters for frame-based processing
    window_dur        = 16; % duration of overlapping frame in second
    step_size_dur     = 1; % frame step-size usually 1 second
    FFT_res_factor    = 2000; % FFT resolution = 1/FFT_res_factor Hz
    f_ref{i}          = func_STFT_single_tone(ref,fs_ref,window_dur,step_size_dur,50,0.1,FFT_res_factor)*2; % reference
    %% 1. single-tone estimation (2nd harmonic)
    f_single{i}       = func_STFT_single_tone(input,FS,window_dur,step_size_dur,fc(1),bound(1),FFT_res_factor);
    %% 2. single-tone enhancement (2nd harmonic)
    f_E_single{i}     = func_STFT_single_tone(input_E_single,FS,window_dur,step_size_dur,fc(1),bound(1),FFT_res_factor);
    %% 3. implementation of [1]: search within sum of harmonic components, mapped to 2nd harmonic
    % for fair comparison use doubled FFT length to ensure consistent search grid at 2nd harmonic.
    f_MLE{i}          = func_STFT_multi_tone_search(input,FS,window_dur,step_size_dur,fc,bound,2*FFT_res_factor);
    %% 4. implementation of [2]: search within weighted sum of harmonic components, mapped to 2nd harmonic
    % for fair comparison use doubled FFT length to ensure consistent search grid at 2nd harmonic.
    f_WMLE{i}         = func_STFT_multi_tone_search_weighted(input,FS,window_dur,step_size_dur,fc,bound,2*FFT_res_factor);
    %% 5. only enhancement MLE
    f_E_MLE{i}        = func_STFT_multi_tone_search(input_E_multi,FS,window_dur,step_size_dur,fc,bound,2*FFT_res_factor);
    %% 6. only enhancement WMLE
    f_E_WMLE{i}       = func_STFT_multi_tone_search_weighted(input_E_multi,FS,window_dur,step_size_dur,fc,bound,2*FFT_res_factor);
    %% 7. only harmonic selection MLE
    correlation_thre    = func_empirical_coeff_thre(length(f_ref{i}),10000)*4; % empirical correlation threshold
    correlation_thre    = min(correlation_thre, 0.8);
    [f_S_MLE{i}, ~, selected_harmonic_index0{i}] = func_STFT_multi_tone_MWC(...
        input,FS,window_dur,step_size_dur,fc,bound,FFT_res_factor,correlation_thre,0);
    %% 8. only harmonic selection MLE
    [f_S_WMLE{i}, ~, selected_harmonic_index0{i}] = func_STFT_multi_tone_MWC(...
        input,FS,window_dur,step_size_dur,fc,bound,FFT_res_factor,correlation_thre,1);
    %% 4. proposed harmonic selection algorithm
    [f_P_MLE{i}, I, selected_harmonic_index{i}] = func_STFT_multi_tone_MWC(...
        input_E_multi,FS,window_dur,step_size_dur,fc,bound,FFT_res_factor,correlation_thre,0);
    f_P_WMLE{i}                                  = func_STFT_multi_tone_MWC(...
        input_E_multi,FS,window_dur,step_size_dur,fc,bound,FFT_res_factor,correlation_thre,1);
    MSE_single(i)       = 1/length(f_ref{i})*norm(f_single{i}-f_ref{i}).^2;
    MSE_E_single(i)     = 1/length(f_ref{i})*norm(f_E_single{i}-f_ref{i}).^2;
    MSE_MLE(i)          = 1/length(f_ref{i})*norm(f_MLE{i}-f_ref{i}).^2;
    MSE_WMLE(i)         = 1/length(f_ref{i})*norm(f_WMLE{i}-f_ref{i}).^2;
    MSE_E_MLE(i)        = 1/length(f_ref{i})*norm(f_E_MLE{i}-f_ref{i}).^2;
    MSE_E_WMLE(i)       = 1/length(f_ref{i})*norm(f_E_WMLE{i}-f_ref{i}).^2;
    MSE_S_MLE(i)        = 1/length(f_ref{i})*norm(f_S_MLE{i}-f_ref{i}).^2;
    MSE_S_WMLE(i)       = 1/length(f_ref{i})*norm(f_S_WMLE{i}-f_ref{i}).^2;
    MSE_P_MLE(i)        = 1/length(f_ref{i})*norm(f_P_MLE{i}-f_ref{i}).^2;
    MSE_P_WMLE(i)       = 1/length(f_ref{i})*norm(f_P_WMLE{i}-f_ref{i}).^2;
end
toc;

[~, order] = sort(MSE_P_MLE,'descend');
figure(111);plot(1:130, MSE_single(order),...
    1:130, MSE_E_single(order),...
    1:130, MSE_MLE(order),...
    1:130, MSE_WMLE(order),...
    1:130, MSE_P_MLE(order));
set(gca, 'YScale', 'log');
grid on; axis([0 131 10^(-6) 10^(-1)]);
xlabel('Recording Index');ylabel('MSE');
leg = legend('$f_{\rm{single}}[l]$','$f_{\rm{E-single}}[l]$',...
             '$f_{\rm{MLE}}[l]$','$f_{\rm{WMLE}}[l]$','$f_{\rm{P-MLE}}[l]$');
set(leg,'Interpreter','latex');


[~, order] = sort(MSE_P_MLE,'descend');
figure(111);plot(1:130, MSE_single(order),...
    1:130, MSE_E_single(order),...
    1:130, MSE_P_MLE(order),...
    1:130, MSE_AMTC(order)+eps);
set(gca, 'YScale', 'log');
grid on; axis([0 131 10^(-6) 3*10^(-2)]);
xlabel('Recording Index');ylabel('MSE');
leg = legend('$f_{\rm{single}}[l]$','$f_{\rm{E-single}}[l]$','$f_{\rm{P-MLE}}[l]$','AMTC');
set(leg,'Interpreter','latex');

figure(112);plot(1:130, MSE_MLE(order),...
    1:130, MSE_WMLE(order),...
    1:130, MSE_P_MLE(order));
set(gca, 'YScale', 'log');
grid on; axis([0 131 10^(-6) 3*10^(-2)]);
xlabel('Recording Index');ylabel('MSE');
leg = legend('$f_{\rm{MLE}}[l]$','$f_{\rm{WMLE}}[l]$','$f_{\rm{P-MLE}}[l]$');
set(leg,'Interpreter','latex');


% figure(1);plot(1:130, MSE_multi, 1:130, MSE_MWC_MLE,1:130, MSE_single,'k--');
% set(gca, 'YScale', 'log');

mean(MSE_single)
mean(MSE_E_single)
mean(MSE_MLE)
mean(MSE_WMLE)
mean(MSE_E_MLE)
mean(MSE_E_WMLE)
mean(MSE_S_MLE)
mean(MSE_S_WMLE)
mean(MSE_P_MLE)
mean(MSE_P_WMLE)

std(MSE_single)
std(MSE_E_single)
std(MSE_MLE)
std(MSE_WMLE)
std(MSE_E_MLE)
std(MSE_E_WMLE)
std(MSE_S_MLE)
std(MSE_S_WMLE)
std(MSE_P_MLE)
std(MSE_P_WMLE)

% 
% 
% nn=zeros(1,130);
% for i = 1:130
%     nn(i) = length(selected_harmonic_index0{i});
% end



% 
% MSE_AMTC=zeros(1,130);
% for i = 1:130
%     e_rref = enf_ref_doc{i}';
%     e_eest = enf_est_doc{i};
%     MSE_AMTC(i) = (sum((e_rref-e_eest).^2))/(length(e_rref));
% end
% figure;plot(1:130,mse_doc,1:130,MSE_AMTC)



% for i = 1:130
%     O1(i) = length(selected_harmonic_index{i});
%     O2(i) = length(selected_harmonic_index0{i});
% end
% [~, s_index] = sort(O2,'descend');
% figure;plot(1:130,O2(s_index),1:130,O1(s_index));
% axis([0 132 0 7]);grid on;
% lll = legend('Without HRFA','With HRFA');
% xxx = xlabel('ENF-WHU Sorted Sample Index');
% yyy = ylabel('$|\Omega|$');
% set(xxx,'Interpreter','latex');
% set(yyy,'Interpreter','latex');
% set(lll,'Interpreter','latex');

















