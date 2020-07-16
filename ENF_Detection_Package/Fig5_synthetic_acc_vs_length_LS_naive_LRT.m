%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program calculates ENF detection accuracies versus observation
% length for the least squares (LS)-likelihood ratio test (LRT) and 
% naive-LRT using synthetic data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
fs        = 400;
T         = 1/fs;
L         = 1000;    % Set to 10000 to generate more similar results as Fig.5
SNR       = -25;                          
duration  = 5:5:120;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load empirically obtained threshold info %%%
load Threshold_info_-25dB_5-5-200_TS2_TS3_10000;
% mean20 : mean of test statistics under H0 for LS-LRT
% mean21 : mean of test statistics under H1 for LS-LRT
% mean30 : mean of test statistics under H0 for naive-LRT
% mean31 : mean of test statistics under H1 for naive-LRT
% var20 : var of test statistics under H0 for LS-LRT
% var21 : var of test statistics under H1 for LS-LRT
% var30 : var of test statistics under H0 for naive-LRT
% var31 : var of test statistics under H1 for naive-LRT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bandpass Filter
F = [0  0.2 0.248 0.249  0.25  0.251 0.252 0.5 1];
M = [0  0   0     0.2      1     0.2     0     0   0 ];
BPF = fir2(255,F,M);
BPFF = abs(fft(BPF,8192));
scalar = max(BPFF);
BPF = BPF/scalar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ACC20          = zeros(1,length(duration));  % accuracy of LS-LRT alpha = 0
ACC21          = zeros(1,length(duration));  % accuracy of LS-LRT alpha = 1
ACC22          = zeros(1,length(duration));  % accuracy of LS-LRT alpha = 2
ACC30          = zeros(1,length(duration));  % accuracy of naive-LRT alpha = 0
ACC31          = zeros(1,length(duration));  % accuracy of naive-LRT alpha = 1
ACC32          = zeros(1,length(duration));  % accuracy of naive-LRT alpha = 2

for j = 1:length(duration)
    j
    result20         = zeros(1,L);
    result21         = zeros(1,L);
    result22         = zeros(1,L);
    result30         = zeros(1,L);
    result31         = zeros(1,L);
    result32         = zeros(1,L);
    Test_Statistic2 = zeros(1,L);
    Test_Statistic3 = zeros(1,L);
    
    ground_truth    = randi(2,[1,L])-1;
    
    duration1  = duration(j);
    N         = fs*duration1;
    ENF      = zeros(1,N);
    A        = 1 + randn(1,N)*0.005;
    while(1)
        f0       = randn(1,N);
        f        = filter(1,[1,-1],f0)*0.0005 + 50;
        if var(f) >= 4*10^(-4) && var(f) <= 5*10^(-4)  % Select appropriate synthetic ENF 
            break;
        end
    end
    phi    = random('unif',0,2*pi,1,1);
    for n = 1:N
        ENF(n)  = A(n)*cos(2*pi/fs*sum(f(1:n)) + phi);
    end
    ENF     =  ENF / norm(ENF);
    
    for count = 1:L
        w        = randn(1,N);
        w        = w / norm(w);
        w        = w ./ (10^(SNR/20));
        
        if ground_truth(count) == 0
            x = w;
        end
        if ground_truth(count) == 1
            x  = ENF + w;
        end
        x_filtered    = filter(BPF,1,x);
        
        
        NFFT_full    = max(2^18,2^(nextpow2(N)+2));
        X_filtered   = abs(fft(x_filtered,NFFT_full));
        X_filtered   = X_filtered(1:(end/2+1));
        
        fc           = find(X_filtered==max(X_filtered))*(fs/NFFT_full);
        Hc           = [cos(2*pi*T*fc*(0:N-1))',sin(2*pi*T*fc*(0:N-1))'];
        Test_Statistic2(count) = 2/N*(x_filtered*Hc)*(Hc'*x_filtered')/((norm(x_filtered).^2));
        
        H_naive      = [cos(2*pi*T*50*(0:N-1))',sin(2*pi*T*50*(0:N-1))'];
        Test_Statistic3(count) = 2/N*(x_filtered*H_naive)*(H_naive'*x_filtered')/((norm(x_filtered).^2));
        
        
        if Test_Statistic2(count) >= mean20(j)+0*sqrt(var20(j))
            result20(count) = 1;
        end
        if Test_Statistic2(count) >= mean20(j)+1*sqrt(var20(j))
            result21(count) = 1;
        end
        if Test_Statistic2(count) >= mean20(j)+2*sqrt(var20(j))
            result22(count) = 1;
        end
        
        if Test_Statistic3(count) >= mean30(j)+0*sqrt(var30(j))
            result30(count) = 1;
        end
        if Test_Statistic3(count) >= mean30(j)+1*sqrt(var30(j))
            result31(count) = 1;
        end
        if Test_Statistic3(count) >= mean30(j)+2*sqrt(var30(j))
            result32(count) = 1;
        end
        
        
    end
    
    [O_TP,O_TN,O_FP,O_FN] = fun_TP_TN_FP_FN(result20,ground_truth);
    ACC20(j)               = (O_TP+O_TN)/L;
    [O_TP,O_TN,O_FP,O_FN] = fun_TP_TN_FP_FN(result21,ground_truth);
    ACC21(j)               = (O_TP+O_TN)/L;
    [O_TP,O_TN,O_FP,O_FN] = fun_TP_TN_FP_FN(result22,ground_truth);
    ACC22(j)    = (O_TP+O_TN)/L;
    
    [O_TP,O_TN,O_FP,O_FN] = fun_TP_TN_FP_FN(result30,ground_truth);
    ACC30(j)               = (O_TP+O_TN)/L;
    [O_TP,O_TN,O_FP,O_FN] = fun_TP_TN_FP_FN(result31,ground_truth);
    ACC31(j)               = (O_TP+O_TN)/L;
    [O_TP,O_TN,O_FP,O_FN] = fun_TP_TN_FP_FN(result32,ground_truth);
    ACC32(j)               = (O_TP+O_TN)/L;
    
end
figure(1);grid on;
plot(duration,ACC20,'ro-',duration,ACC21,'rx-',duration,ACC22,'r+-',...
    duration,ACC30,'ko-',duration,ACC31,'kx-',duration,ACC32,'k+-');
grid on; axis([0 125 0 1.02]);
hl = legend('LS-LRT, $\alpha = 0$','LS-LRT, $\alpha = 1$','LS-LRT, $\alpha = 2$',...
            'naive-LRT, $\alpha = 0$','naive-LRT, $\alpha = 1$','naive-LRT, $\alpha = 2$');
hx = xlabel('$N/f_{\rm{S}}$');
hy = ylabel('Accuracy ($\%$)');
set(hx, 'Interpreter', 'latex');
set(hy, 'Interpreter', 'latex');
set(hl, 'Interpreter', 'latex');
% T00 = Test_Statistic3((ground_truth==0));
% T01 = Test_Statistic3((ground_truth==1));
% figure(12);histogram(T00,100);hold on; histogram(T01,100);






















