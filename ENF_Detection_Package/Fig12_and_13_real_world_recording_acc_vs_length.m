%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program calculates ENF detection accuracies of 
%
% 2. LS-LRT,
% 3. naive-LRT, 
% 4. TF detector,
%
% versus recording length using real-world audio recordings.
%
% Index "1" is used for matched filter which is unavailable here.
%
%   Note: the results presented in the paper used 50 recordings (by the time
%   of the preparation of the paper) under H1 and 10 under H0 (extended to 50
%   via random cropping). The current full dataset contains 130 recordings 
%   under H1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bandpass Filter
%%% length of BPF increased to 1024 to tackle real-world recordings
F2 = [0 0.4 0.499 0.4995 0.5 0.5005 0.501 0.6 0.8 1];
M2 = [0 0 0 0.2 1 0.2 0 0 0 0];
BPF= fir2(1023,F2,M2);
BPFF     = abs(fft(BPF,8192));
scalar   = max(BPFF);
BPF      = BPF/scalar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs              = 400;
T               = 1/fs;
AWindowLength   = 16*fs;
AWindowShift    = rectwin(AWindowLength)';
AStepSize       = 1*fs;
NFFT            = 200*fs;

%%%%% setting for Fig. 12 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
duration = 10:20:270;
load Threshold_info_10_20_270
thre2 = mean20+2*sqrt(var20);
thre3 = mean30+2*sqrt(var30);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% setting for Fig. 13 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% duration = 5:10;
% load Threshold_info_5_1_10
% thre2 = mean20+2*sqrt(var20);
% thre3 = mean30+2*sqrt(var30);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% change path to where you store the recordings
path         = 'Here is your path of recordings';
H0_index     = dir(strcat(path,'H0'));
H1_index     = dir(strcat(path,'H1'));
ground_truth = [ones(1,length(H1_index)-2),zeros(1,length(H0_index)-2)];

% result1 omitted since it is to denote matched filter, not applicable here
result2      = zeros(length(duration),(length(H1_index)-2+length(H0_index)-2));   % LS-LRT
result3      = zeros(length(duration),(length(H1_index)-2+length(H0_index)-2));   % naive-LRT
result4      = zeros(length(duration),(length(H1_index)-2+length(H0_index)-2));   % TF
ACC2         = zeros(1,length(duration));
ACC3         = zeros(1,length(duration));
ACC4         = zeros(1,length(duration));

O_TP2 = zeros(1,length(duration));
O_TN2 = zeros(1,length(duration));
O_FP2 = zeros(1,length(duration));
O_FN2 = zeros(1,length(duration));
O_TP3 = zeros(1,length(duration));
O_TN3 = zeros(1,length(duration));
O_FP3 = zeros(1,length(duration));
O_FN3 = zeros(1,length(duration));
O_TP4 = zeros(1,length(duration));
O_TN4 = zeros(1,length(duration));
O_FP4 = zeros(1,length(duration));
O_FN4 = zeros(1,length(duration));

for i = 1:(length(H1_index)-2+length(H0_index)-2)
    i
    if i<=(length(H1_index)-2)
        [audio, fs0] = audioread(strcat(H1_index(i+2).folder,'\',H1_index(i+2).name));
        audio        = audio(:,1)';
    else
        [audio, fs0] = audioread(strcat(H0_index((i-(length(H1_index)-2))+2).folder,'\',H0_index((i-(length(H1_index)-2))+2).name));
        audio        = audio(:,1)';
    end
    
    for j = 1:length(duration)
        current_dur  = duration(j);
        start_index  = randi(length(audio)-current_dur*fs0);
        audio_cut    = audio(start_index:(start_index+current_dur*fs0-1));
        x            = resample(audio_cut, fs, fs0);
        N            = length(x);
        
        x_filtered   = filter(BPF,1,x);
        
        NFFT_full    = max(2^18,2^(nextpow2(N)+2));
        X_filtered   = abs(fft(x_filtered,NFFT_full));
        X_filtered   = X_filtered(1:(end/2+1));
        
        fc              = find(X_filtered==max(X_filtered))*(fs/NFFT_full);
        Hc              = [cos(2*pi*T*fc*(0:N-1))',sin(2*pi*T*fc*(0:N-1))'];
        Test_Statistic2 = 2/N*(x_filtered*Hc)*(Hc'*x_filtered')/((norm(x_filtered).^2));
        
        H_naive         = [cos(2*pi*T*100*(0:N-1))',sin(2*pi*T*100*(0:N-1))'];
        Test_Statistic3 = 2/N*(x_filtered*H_naive)*(H_naive'*x_filtered')/((norm(x_filtered).^2));
        
        IF              = fun_STFT_interpo(x_filtered, AWindowShift, AStepSize,fs,NFFT);
        Test_Statistic4 = var(IF);
        
        
        if Test_Statistic2 >= thre2(j)
            result2(j,i) = 1;
        end
        if Test_Statistic3 >= thre3(j)
            result3(j,i) = 1;
        end
        if Test_Statistic4 < 0.005 % threshold according to new BPF 
            result4(j,i) = 1;
        end
    end
end
for j = 1:length(duration)
    [O_TP2(j),O_TN2(j),O_FP2(j),O_FN2(j)] = fun_TP_TN_FP_FN(result2(j,:),ground_truth);
    ACC2(j)                  = (O_TP2(j)+O_TN2(j))/(length(H1_index)-2+length(H0_index)-2);
    [O_TP3(j),O_TN3(j),O_FP3(j),O_FN3(j)] = fun_TP_TN_FP_FN(result3(j,:),ground_truth);
    ACC3(j)                  = (O_TP3(j)+O_TN3(j))/(length(H1_index)-2+length(H0_index)-2);
    [O_TP4(j),O_TN4(j),O_FP4(j),O_FN4(j)] = fun_TP_TN_FP_FN(result4(j,:),ground_truth);
    ACC4(j)                  = (O_TP4(j)+O_TN4(j))/(length(H1_index)-2+length(H0_index)-2);
end

%%%%%%% plot for Fig.%%%%%%% 12%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
plot(duration,ACC2*100,'bo-',duration,ACC3*100,'bx-',duration,ACC4*100,'ro-');
grid on;axis([0 280 0 100]);
hl = legend('LS-LRT','naive-LRT','TF');
hx = xlabel('$N/f_{\rm{S}}$');
hy = ylabel('Accuracy ($\%$)');
set(hx, 'Interpreter', 'latex');
set(hy, 'Interpreter', 'latex');
set(hl, 'Interpreter', 'latex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% plot for Fig.%%%%%%% 13%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2);
% plot(duration,ACC2*100,'bo-',duration,ACC3*100,'bx-',duration,ACC4*100,'ro-');
% grid on;axis([4.8 10.2 0 100]);
% hl = legend('LS-LRT','naive-LRT','TF');
% hx = xlabel('$N/f_{\rm{S}}$');
% hy = ylabel('Accuracy ($\%$)');
% set(hx, 'Interpreter', 'latex');
% set(hy, 'Interpreter', 'latex');
% set(hl, 'Interpreter', 'latex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







