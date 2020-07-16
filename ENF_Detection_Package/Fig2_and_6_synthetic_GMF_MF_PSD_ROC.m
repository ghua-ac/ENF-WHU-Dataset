%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program generates the ROC curves of ENF detectors, including
% GMF, MF approximation, and the asymptotic approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
fs        = 400;                  % Hz
duration  = 5;                    % second
N         = fs*duration;          % sample length
L         = 100;                  % number of random realizations (set to 1000 for paper, here set to 100 to save time)
SNR       = -25;                  % in dB

ENF      = zeros(1,N);
f0       = randn(1,N);
A        = 1 + randn(1,N)*0.005;                  % amplitude
f        = filter(1,[1,-0.99],f0)*0.0005 + 50;    % instantaneous frequency
theta    = random('unif',0,2*pi,1,1);             % phase
for n = 1:N
    ENF(n)  = A(n)*cos(2*pi/fs*sum(f(1:n)) + theta); % clean ENF
end
ENF     =  ENF / norm(ENF);                        % normalization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bandpass Filter
F = [0  0.2 0.248 0.249  0.25  0.251 0.252 0.5 1];
M = [0  0   0     0.2      1     0.2     0     0   0 ];
BPF = fir2(255,F,M);
BPFF = abs(fft(BPF,8192));
scalar = max(BPFF);
BPF = BPF/scalar;

%Show Fig. 6 in paper%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
FR = 20*log10(abs(fft(BPF,4000)));
figure(999); plot(0:0.1:199.9,FR(1:2000));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma           = [-2, -0.5, 0, 0.25, 0.5, 0.75, 1, 1.5, 2.5]; % threshold values acorrding to mean 0 and 1

FPR1          = zeros(1,length(gamma));
TPR1          = zeros(1,length(gamma));
FPR2          = zeros(1,length(gamma));
TPR2          = zeros(1,length(gamma));
FPR3          = zeros(1,length(gamma));
TPR3          = zeros(1,length(gamma));

for j = 1:length(gamma)
    gamma1          = gamma(j)
    result1         = zeros(1,L);        % GMF
    result2         = zeros(1,L);        % MF
    result3         = zeros(1,L);        % Asymptotic
    Test_Statistic1 = zeros(1,L);
    Test_Statistic2 = zeros(1,L);
    Test_Statistic3 = zeros(1,L);
    ground_truth    = randi(2,[1,L])-1;
    for count = 1:L
       count
        w        = randn(1,N);               % noise initialization
        w        = w / norm(w);              % normalization
        w        = w ./ (10^(SNR/20));       % scale noise according to SNR
        
        if ground_truth(count) == 0
            x = w;
        end
        if ground_truth(count) == 1
            x  = ENF + w;
        end
        x_filtered    = filter(BPF,1,x);
        ENF_filtered  = filter(BPF,1,ENF);
        w_filtered    = filter(BPF,1,w);
        
        cov_coeff_real = 1/N*xcorr(w_filtered);
        C = zeros(N,N);
        C_orth = zeros(N,N);
        for ii = 0:N-1
            C(ii+1,:) = cov_coeff_real(((end+1)/2-ii):((end+1)/2-ii)+(N-1));
        end
        Cinv = C^(-1);
         
        %%%%%%%% All test statistics normalized to mean 0 and mean 1 distributions
        
        Test_Statistic1(count) =  (x_filtered*Cinv*ENF_filtered')/(ENF_filtered*Cinv*ENF_filtered'); % GMF
        Test_Statistic2(count) =  (x_filtered*ENF_filtered')/(ENF_filtered*ENF_filtered');  % MF

        nfft          = 2^nextpow2(length(x));
        BPF_freq      = fft(BPF,nfft)/sqrt(nfft);
        ENF_freq      = fft(ENF,nfft)/sqrt(nfft);
        X             = fft(x,nfft)/sqrt(nfft);

        X_Filtered    = fft(x_filtered,nfft)/sqrt(nfft);
        W_Filtered    = fft(w_filtered,nfft)/sqrt(nfft);
        ENF_Filtered  = fft(ENF_filtered,nfft)/sqrt(nfft);

        TS3                    = real(sum((X_Filtered  .*conj(ENF_Filtered))./(W_Filtered.*conj(W_Filtered)))/nfft);
        Test_Statistic3(count) = TS3/(real(sum((ENF_Filtered.*conj(ENF_Filtered))./(W_Filtered.*conj(W_Filtered))))/nfft);  % Asymptotic
        
        %
        if Test_Statistic1(count) >= gamma1
            result1(count) = 1;
        end
        if Test_Statistic2(count) >= gamma1
            result2(count) = 1;
        end
        if Test_Statistic3(count) >= gamma1
            result3(count) = 1;
        end

    end
    [O_TP,O_TN,O_FP,O_FN] = fun_TP_TN_FP_FN(result1,ground_truth);
    FPR1(j)  =   O_FP/(O_FP+O_TN);
    TPR1(j)  =   O_TP/(O_TP+O_FN);
    [O_TP,O_TN,O_FP,O_FN] = fun_TP_TN_FP_FN(result2,ground_truth);
    FPR2(j)  =   O_FP/(O_FP+O_TN);
    TPR2(j)  =   O_TP/(O_TP+O_FN);
    [O_TP,O_TN,O_FP,O_FN] = fun_TP_TN_FP_FN(result3,ground_truth);
    FPR3(j)  =   O_FP/(O_FP+O_TN);
    TPR3(j)  =   O_TP/(O_TP+O_FN);
end
figure(1);hold on;grid on;
plot(0:0.1:1,0:0.1:1,'--',FPR1,TPR1,'ro-',FPR2,TPR2,'bo-',FPR3,TPR3,'ko-');
hl = legend('','GMF (7)', 'MF (9)','Asymptotic (17)');
hx = xlabel('$P_{\rm{FA}}$');
hy = ylabel('$P_{\rm D}$');
set(hx, 'Interpreter', 'latex');
set(hy, 'Interpreter', 'latex');
set(hl, 'Interpreter', 'latex');

T00 = Test_Statistic3((ground_truth==0));
T01 = Test_Statistic3((ground_truth==1));
figure(11);histogram(T00,100);hold on; histogram(T01,100);
% % var(T00)
% % var(T01)



