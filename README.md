# About

This GitHub link contains the ENF-WHU audio recording dataset collected around Wuhan University campus, and the MATLAB programs for electronic network frequency (ENF) analysis in ENF-based audio forensic applications.

## ENF-WHU Dataset
- **Recording location:** classroom, campus path, meeting room, graduate student office, dormitory, library.
- **Environment diversity:** day/night, rainy/suny.
- **Recording device:** popular smartphone and voice recorder.
- **Duration:** 5~20 minutes
- **Format:** PCM WAVE
- **Quantization depth:** 16-bit
- **Channel:** mono
- **Sampling frequencuy:** 8000 Hz
- **Category:**
  1 **H0:** "O01\~O10.wav" 10 real-world recordings without captured ENF. "01\~40.wav" 40 segments under H0 are obtained by random cropping the 10 recordings. 
  2 **H1:** "001~130.wav" 130 real-world recordings with captured (noisy) ENF.
  3 **H1_ref:** "001_ref~130_ref.wav" the corresponding 130 reference ENF (noise-free) obtained from power main.

## Citation Information
- ENF Detection:
 > \[1] G. Hua, H. Liao, Q. Wang, H. Zhang, and D. Ye, "Detection of electric network frequency in audio recordings – From theory to practical detectors," IEEE Trans. Inf. Forensics Security, DOI: 10.1109/TIFS.2020.3009579, 2020.
- ENF Enhancement:
 > \[2] G. Hua and H. Zhang, “ENF signal enhancement in audio recordings,” IEEE Transactions on Information Forensics and Security, vol. 15, pp. 1868-1878, 2020.

## MATLAB Programs
- ENF Detection
- ENF Enhancement
