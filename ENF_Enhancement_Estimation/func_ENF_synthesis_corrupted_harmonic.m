function [raw_wave, IF_2nd, raw_wave_multi] = func_ENF_synthesis_corrupted_harmonic(fundamental_f, harmonic_index, corrupted_index, duration, fs)
N               = duration*fs;
%% create fundamental IF
f0              = randn(1,N);
IF0             = filter(1,[1,-1],f0)*0.0005;
IF0             = IF0/std(IF0)*sqrt(4.5*10^(-4));
IF0             = IF0 + fundamental_f;
IFs             = harmonic_index'*IF0; % IFs across all harmonics
index           = intersect(harmonic_index,corrupted_index)-1;
for i = 1:length(index)
      IFs(index(i),:) = IFs(index(i),:) + 5*randn(1,N);
end
%% instantaneous amplutudes and initial phases
% figure(998);hold on; plot(IFs(1,:));plot(IFs(2,:));plot(IFs(3,:));plot(IFs(4,:));plot(IFs(5,:));plot(IFs(6,:));
% figure(999);hold on; plot(IFs(1,:));
N_harmonic      = length(harmonic_index);
amps            = 1 + randn(N_harmonic,N)*0.005; % instantaneous amplitudes
phases          = random('unif',0,2*pi,N_harmonic,1); % initial phases
%% synthesize time domain waveforms
ENF_multi = zeros(N_harmonic, N);
for n = 1:N
    ENF_multi(:,n) = amps(:,n).*cos(2*pi/fs*sum(IFs(:,1:n),2) + phases); % this is very time-consuming
end
for i = 1:6
    ENF_multi(i,:) = ENF_multi(i,:)/norm(ENF_multi(i,:));
end
raw_wave         = sum(ENF_multi);
raw_wave_multi   = ENF_multi;
raw_wave         = raw_wave / norm(raw_wave); % ensure unit norm
IF_2nd           = 2*IF0;
