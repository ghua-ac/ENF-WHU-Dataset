function [IF0,IF1] =  fun_STFT_interpo(Signal,Window,StepPoints,Fs,NFFT)

FrameSize       = length(Window);
% Signal = [Signal  zeros(1,Fs - mod(length(Signal),Fs))];
Signal          = [zeros(1, FrameSize/2), Signal, zeros(1, FrameSize/2)];
WindowPositions = 1: StepPoints: length(Signal) - FrameSize+1;
for i = 1:length(WindowPositions)
   temp = fft(Signal(WindowPositions(i):WindowPositions(i) + FrameSize - 1).*Window,NFFT);
   HalfTempFFT  = temp(1:end/2);
   absHalfTempFFT = abs(HalfTempFFT);
   ValueMax     = max(absHalfTempFFT);
   PeakLoc      = find(absHalfTempFFT==ValueMax(1));
   ValueLeft    = HalfTempFFT(PeakLoc - 1);
   ValueCenter  = HalfTempFFT(PeakLoc);
   ValueRight   = HalfTempFFT(PeakLoc + 1);
   CorrectionCoef = -real((ValueRight - ValueLeft)/(2*ValueCenter - ValueRight - ValueLeft));
   IF0(i)       = (PeakLoc + CorrectionCoef -1)*Fs/NFFT;
   IF1(i)       = (PeakLoc - 1)*Fs/NFFT;
end
