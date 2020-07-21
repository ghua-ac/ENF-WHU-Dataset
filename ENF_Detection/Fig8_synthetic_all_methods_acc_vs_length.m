%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program calculates ENF detection accuracies versus observation
% length for
%
% 1. MF
% 2. LS-LRT,
% 3. naive-LRT, 
% 4. TF detector,
%
% using synthetic data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;
fs        = 400;
T         = 1/fs;                     
L         = 100;    % Set to 1000 to generate results similar to Fig. 8
SNR       = -25;                        
duration  = 5:5:200;
load Threshold_info_-25dB_5-5-200_TS2_TS3_10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bandpass Filter
F = [0  0.2 0.248 0.249  0.25  0.251 0.252 0.5 1];
M = [0  0   0     0.2      1     0.2     0     0   0 ];
BPF = fir2(255,F,M);
BPFF = abs(fft(BPF,8192));
scalar = max(BPFF);
BPF = BPF/scalar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ACC1          = zeros(1,length(duration));
ACC2          = zeros(1,length(duration));
ACC3          = zeros(1,length(duration));
ACC4          = zeros(1,length(duration));

for j = 1:length(duration)
    j
    result1         = zeros(1,L);
    result2         = zeros(1,L);
    result3         = zeros(1,L);
    result4         = zeros(1,L);
    Test_Statistic1 = zeros(1,L);
    Test_Statistic2 = zeros(1,L);
    Test_Statistic3 = zeros(1,L);
    Test_Statistic4 = zeros(1,L);
    ground_truth    = randi(2,[1,L])-1;
    
    duration1  = duration(j);
    N         = fs*duration1;
    ENF      = zeros(1,N);
    A        = 1 + randn(1,N)*0.005;
    while(1)
        f0       = randn(1,N);
        f        = filter(1,[1,-1],f0)*0.0005 + 50;
        if var(f) >= 4*10^(-4) && var(f) <= 5*10^(-4)
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
        ENF_filtered  = filter(BPF,1,ENF);
        
        Test_Statistic1(count) =  (x_filtered*ENF_filtered');
        
        NFFT_full    = max(2^18,2^(nextpow2(N)+2));
        X_filtered   = abs(fft(x_filtered,NFFT_full));
        X_filtered   = X_filtered(1:(end/2+1));

        fc           = find(X_filtered==max(X_filtered))*(fs/NFFT_full);
        Hc           = [cos(2*pi*T*fc*(0:N-1))',sin(2*pi*T*fc*(0:N-1))']; 
        Test_Statistic2(count) = 2/N*(x_filtered*Hc)*(Hc'*x_filtered')/((norm(x_filtered).^2));
        
        H_naive      = [cos(2*pi*T*50*(0:N-1))',sin(2*pi*T*50*(0:N-1))'];
        Test_Statistic3(count) = 2/N*(x_filtered*H_naive)*(H_naive'*x_filtered')/((norm(x_filtered).^2));

        
        AWindowLength1 = 16*fs;
        AWindowShift1 = rectwin(AWindowLength1)';
        AStepSize = 1*fs;
        NFFT = 200*fs;
        
        IF=fun_STFT_interpo(x_filtered, AWindowShift1, AStepSize,fs,NFFT);
        Test_Statistic4(count) = var(IF);          
                
        if Test_Statistic1(count) >= 0.5  % according to mean 0 and mean 1
            result1(count) = 1;
        end
        if Test_Statistic2(count) >= mean20(j)+2*sqrt(var20(j))
            result2(count) = 1;
        end
        if Test_Statistic3(count) >= mean30(j)+2*sqrt(var30(j))
            result3(count) = 1;
        end     
        if Test_Statistic4(count) < 0.08  % empirical value according to mean40 and mean41
            result4(count) = 1;
        end 
    end
    [O_TP,O_TN,O_FP,O_FN] = fun_TP_TN_FP_FN(result1,ground_truth);
    ACC1(j)               = (O_TP+O_TN)/L;
    [O_TP,O_TN,O_FP,O_FN] = fun_TP_TN_FP_FN(result2,ground_truth);
    ACC2(j)               = (O_TP+O_TN)/L;
    [O_TP,O_TN,O_FP,O_FN] = fun_TP_TN_FP_FN(result3,ground_truth);
    ACC3(j)               = (O_TP+O_TN)/L;
    [O_TP,O_TN,O_FP,O_FN] = fun_TP_TN_FP_FN(result4,ground_truth);
    ACC4(j)               = (O_TP+O_TN)/L;
end
figure(1);hold on;grid on;
plot(duration,ACC1*100,'ro-',duration,ACC2*100,'bo-',duration,ACC3*100,'bx-',duration,ACC4*100,'k-');
axis([0 200 50 102]);
hl = legend('MF','LS-LRT','naive-LRT','TF');
hx = xlabel('$N/f_{\rm{S}}$');
hy = ylabel('Accuracy ($\%$)');
set(hx, 'Interpreter', 'latex');
set(hy, 'Interpreter', 'latex');
set(hl, 'Interpreter', 'latex');




