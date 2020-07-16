%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program calculates the ROC curves of ENF detection for
%
% 2. LS-LRT,
% 3. naive-LRT, 
% 4. TF detector,
%
% using using real-world audio recordings from ENF-WHU dataset. Two values 
% of recording lengths are available, i.e., 150 and 5 seconds respectively.
%
% Index "1" is used for matched filter which is unavailable here.
%
%   Note: the results presented in the paper used 50 recordings (by the time
%   of the preparation of the paper) under H1 and 10 under H0 (extended to 50
%   via random cropping). The current full dataset contains 130 recordings 
%   under H1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;

F = [0 0.4 0.499 0.4995 0.5 0.5005 0.501 0.6 0.8 1];
M = [0 0 0 0.2 1 0.2 0 0 0 0];
BPF= fir2(1023,F,M);
BPFF     = abs(fft(BPF,8192));
scalar   = max(BPFF);
BPF      = BPF/scalar;

fs              = 400;
T               = 1/fs;
AWindowLength   = 16*fs;
AWindowShift    = rectwin(AWindowLength)';
AStepSize       = 1*fs;
NFFT            = 200*fs;

% %%%%%%%%%%%%%%%%%%%%%%%%%% threshold values for recording of 150 seconds %%%%%%%%%%%%%%%%%%%%%
% N_thre = 1000;
% duration = 150;
% load Threshold_info_10_20_270
% % thre2 = mean20(10)+[-2*sqrt(var20(10)), 1*sqrt(var20(10)), 3*sqrt(var20(10)),5*sqrt(var20(10)),...
% %     8*sqrt(var20(10)), 15*sqrt(var20(10)), 35*sqrt(var20(10))];
% % thre3 = mean30(10)+[-2*sqrt(var30(10)), 0, 1*sqrt(var30(10)),5*sqrt(var30(10))...
% %     7*sqrt(var30(10)), 10*sqrt(var30(10)), 20*sqrt(var30(10))];
% % thre4 = [0.00005 0.001 0.0015 0.004 0.009 0.014 0.02];
% 
% thre2 = linspace(mean20(10)-2*sqrt(var20(10)),mean20(10)+35*sqrt(var20(10)),N_thre);
% thre3 = linspace(mean30(10)-2*sqrt(var30(10)),mean30(10)+20*sqrt(var30(10)),N_thre);
% thre4 = linspace(0.00005,0.02,N_thre);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% threshold values for recording of 5 seconds %%%%%%%%%%%%%%%%%%%%%
N_thre = 1000;
duration = 5;
load Threshold_info_5_1_10
% thre2 = mean20(1)+[-3*sqrt(var20(1)), -1*sqrt(var20(1)), 0, 1*sqrt(var20(1)),1.5*sqrt(var20(1)),...
%                    2*sqrt(var20(1)), 3*sqrt(var20(1))];
% thre3 = mean30(1)+[-2*sqrt(var30(1)), -1*sqrt(var30(1)), 0,1*sqrt(var30(1))...
%                    1.5*sqrt(var30(1)), 2*sqrt(var30(1)), 3*sqrt(var30(1))];
% thre4 = [10^(-11), 10^(-10), 10^(-9), 10^(-8.5), 10^(-8), 10^(-7), 10^(-6)];

thre2 = linspace(mean20(1)-2*sqrt(var20(1)),mean20(1)+35*sqrt(var20(1)),N_thre);
thre3 = linspace(mean30(1)-2*sqrt(var30(1)),mean30(1)+20*sqrt(var30(1)),N_thre);
thre4 = linspace(10^(-11),10^(-6),N_thre);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path         = 'Here is your path of recordings';
H0_index     = dir(strcat(path,'H0'));
H1_index     = dir(strcat(path,'H1'));
ground_truth = [ones(1,length(H1_index)-2),zeros(1,length(H0_index)-2)];

Test_Statistic2 = zeros(1,length(H1_index)-2+length(H0_index)-2);
Test_Statistic3 = zeros(1,length(H1_index)-2+length(H0_index)-2);
Test_Statistic4 = zeros(1,length(H1_index)-2+length(H0_index)-2);
result2      = zeros((N_thre),(length(H1_index)-2+length(H0_index)-2));
result3      = zeros((N_thre),(length(H1_index)-2+length(H0_index)-2));
result4      = zeros((N_thre),(length(H1_index)-2+length(H0_index)-2));
ACC2         = zeros(1,(N_thre));
ACC3         = zeros(1,(N_thre));
ACC4         = zeros(1,(N_thre));

O_TP2 = zeros(1,(N_thre));
O_TN2 = zeros(1,(N_thre));
O_FP2 = zeros(1,(N_thre));
O_FN2 = zeros(1,(N_thre));
O_TP3 = zeros(1,(N_thre));
O_TN3 = zeros(1,(N_thre));
O_FP3 = zeros(1,(N_thre));
O_FN3 = zeros(1,(N_thre));
O_TP4 = zeros(1,(N_thre));
O_TN4 = zeros(1,(N_thre));
O_FP4 = zeros(1,(N_thre));
O_FN4 = zeros(1,(N_thre));


for i = 1:(length(H1_index)-2+length(H0_index)-2)
    i
    if i<=(length(H1_index)-2)
        [audio, fs0] = audioread(strcat(H1_index(i+2).folder,'\',H1_index(i+2).name));
        audio        = audio(:,1)';
    else
        [audio, fs0] = audioread(strcat(H0_index((i-(length(H1_index)-2))+2).folder,'\',H0_index((i-(length(H1_index)-2))+2).name));
        audio        = audio(:,1)';
    end
    
    current_dur  = duration;
    start_index  = randi(length(audio)-current_dur*fs0);
    audio_cut    = audio(start_index:(start_index+current_dur*fs0-1));
    x            = resample(audio_cut, fs, fs0);
    N            = length(x);
    
    x_filtered   = filter(BPF,1,x);
    
    NFFT_full    = max(2^18,2^(nextpow2(N)+2));
    X_filtered   = abs(fft(x_filtered,NFFT_full));
    X_filtered   = X_filtered(1:(end/2+1));
    
    fc                 = find(X_filtered==max(X_filtered))*(fs/NFFT_full);
    Hc                 = [cos(2*pi*T*fc*(0:N-1))',sin(2*pi*T*fc*(0:N-1))'];
    Test_Statistic2(i) = 2/N*(x_filtered*Hc)*(Hc'*x_filtered')/((norm(x_filtered).^2));
    
    H_naive            = [cos(2*pi*T*100*(0:N-1))',sin(2*pi*T*100*(0:N-1))'];
    Test_Statistic3(i) = 2/N*(x_filtered*H_naive)*(H_naive'*x_filtered')/((norm(x_filtered).^2));
    
    IF                 = fun_STFT_interpo(x_filtered, AWindowShift, AStepSize,fs,NFFT);
    Test_Statistic4(i) = var(IF);
    
    for j = 1:N_thre
        if Test_Statistic2(i) >= thre2(j)
            result2(j,i) = 1;
        end
        if Test_Statistic3(i) >= thre3(j)
            result3(j,i) = 1;
        end
        if Test_Statistic4(i) < thre4(j)
            result4(j,i) = 1;
        end
    end
end

for j = 1:N_thre
    [O_TP2(j),O_TN2(j),O_FP2(j),O_FN2(j)] = fun_TP_TN_FP_FN(result2(j,:),ground_truth);
    ACC2(j)                  = (O_TP2(j)+O_TN2(j))/(length(H1_index)-2+length(H0_index)-2);
    [O_TP3(j),O_TN3(j),O_FP3(j),O_FN3(j)] = fun_TP_TN_FP_FN(result3(j,:),ground_truth);
    ACC3(j)                  = (O_TP3(j)+O_TN3(j))/(length(H1_index)-2+length(H0_index)-2);
    [O_TP4(j),O_TN4(j),O_FP4(j),O_FN4(j)] = fun_TP_TN_FP_FN(result4(j,:),ground_truth);
    ACC4(j)                  = (O_TP4(j)+O_TN4(j))/(length(H1_index)-2+length(H0_index)-2);
end

figure(1);
plot(O_FP2/(length(H0_index)-2),O_TP2/(length(H1_index)-2),'bo-',...
    O_FP3/(length(H0_index)-2),O_TP3/(length(H1_index)-2),'bx-',...
    O_FP4/(length(H0_index)-2),O_TP4/(length(H1_index)-2),'ro-',...
    [0,1],[0,1]);
grid on;
axis([0 1 0 1]);
hl = legend('LS-LRT','naive-LRT','TF');
hx = xlabel('$P_{\rm{FA}}$');
hy = ylabel('$P_{\rm{D}}$');
set(hx, 'Interpreter', 'latex');
set(hy, 'Interpreter', 'latex');
set(hl, 'Interpreter', 'latex');


%% Check the distribution of test statistics
figure(2);
num_histogram = 100;
T21 = Test_Statistic2(1:50);
T20 = Test_Statistic2(51:100);
subplot(311);histogram(T21,num_histogram);hold on; histogram(T20,num_histogram );
legend('H_1','H_0');xlabel('Test Statistic Value');ylabel('Count');title('LS-LRT');
T31 = Test_Statistic3(1:50);
T30 = Test_Statistic3(51:100);
subplot(312);histogram(T31,num_histogram);hold on; histogram(T30,num_histogram);
legend('H_1','H_0');xlabel('Test Statistic Value');ylabel('Count');title('naive-LRT')
T41 = Test_Statistic4(1:50);
T40 = Test_Statistic4(51:100);
subplot(313);histogram(T41,num_histogram);hold on; histogram(T40,num_histogram);
legend('H_1','H_0');xlabel('Test Statistic Value');ylabel('Count');title('TF');

%% calculate AUC using Matlab function
score2 = sum(result2)/N_thre; 
score3 = sum(result3)/N_thre; 
score4 = sum(result4)/N_thre; 
[X2,Y2,T2,AUC2] = perfcurve(ground_truth',score2',1);
[X3,Y3,T3,AUC3] = perfcurve(ground_truth',score3',1);
[X4,Y4,T4,AUC4] = perfcurve(ground_truth',score4',1);

