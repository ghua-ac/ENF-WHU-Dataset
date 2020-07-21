# About

This GitHub link contains the ENF-WHU audio recording dataset collected around Wuhan University campus, and the MATLAB programs for electronic network frequency (ENF) analysis in ENF-based audio forensic applications.

## ENF-WHU Dataset
- **Recording location:** classroom, campus path, meeting room, graduate student office, dormitory, library.
- **Environment diversity:** day/night, rainy/suny, interior/exterior.
- **Recording device:** popular smartphone and voice recorder.
- **Duration:** 5~20 minutes
- **Format:** PCM WAVE
- **Quantization depth:** 16-bit
- **Channel:** mono
- **Sampling frequencuy:** 8000 Hz
- **Category:**<br>
  **H1:** "001\~130.wav" 130 real-world recordings with captured (noisy) ENF.<br>
  **H1_ref:** "001_ref\~130_ref.wav" the corresponding 130 reference ENF (noise-free) obtained from power main.<br>
  **H0:** "O01\~O10.wav" 10 real-world recordings without captured ENF. "01\~40.wav" 40 segments under H0 obtained by random cropping the 10 recordings.

## MATLAB Programs
### ENF Detection
1. Clairvoyant detectors: Neyman-Pearson (NP) detectors assuming perfect knowledge of the ENF.
	1. Generalized matched filter (GMF): a standard NP detector.
	2. Matched filter (MF)-like approximation: avoid the requirement of unknown noise covariance matrix.
	3. Asymptotic approximation: trade-off between computational complexity and detection performance.

2. Generalized likelihood ratio test (GLRT) detectors: ENF assumed unknown and deterministic.
	1. Least-squares (LS)-LRT: MF-like with unknown parameters replaced by the maximum likelihood estimations (MLEs).
	2. naive-LRT: MF-like with the unknown frequencies replaces by the nominal value (50 Hz or harmonics).

3. Time-frequency (TF) domain detector (proposed)
	- Frame-based processing
	- Test statistic is the sample variance of the strongest time-frequency line (e.g., STFT + peak)
	- $$H_0$$


### ENF Enhancement

## Citation Information
- **ENF Detection:**
 > \[1] G. Hua, H. Liao, Q. Wang, H. Zhang, and D. Ye, "Detection of electric network frequency in audio recordings – From theory to practical detectors," IEEE Trans. Inf. Forensics Security, DOI: 10.1109/TIFS.2020.3009579, 2020. [link](https://ieeexplore.ieee.org/document/9143185)
- **ENF Enhancement:**
 > \[2] G. Hua and H. Zhang, "ENF signal enhancement in audio recordings," IEEE Trans. Inf. Forensics Security, vol. 15, pp. 1868-1878, 2020. [link](https://ieeexplore.ieee.org/document/8894138)
- **Related Works:**
 > \[3] G. Hua, "Error analysis of forensic ENF matching," in Proc. 2018 IEEE International Workshop on Information Forensics and Security (WIFS), pp. 1-7, Hong Kong, Dec. 2018. [link](https://ieeexplore.ieee.org/document/8630786)<br>
 > \[4] G. Hua, G. Bi, and V. L. L. Thing, "On practical issues of electric network frequency based audio forensics," IEEE Access, vol. 5, pp. 20640-20651, Oct. 2017. [link](https://ieeexplore.ieee.org/document/7807225)<br>
 > \[5] G. Hua, Y. Zhang, J. Goh, and V. L. L. Thing, "Audio authentication by exploring the absolute error map of the ENF signals," IEEE Trans. Inf. Forensics Security, vol. 11, no. 5, pp. 1003-1016, May 2016. [link](https://ieeexplore.ieee.org/document/7378470)<br>
 > \[6] G. Hua, J. Goh, and V. L. L. Thing, “A dynamic matching algorithm for audio timestamp identification using the ENF criterion,” IEEE Trans. Inf. Forensics Security, vol. 9, no. 7, pp. 1045-1055, Jul. 2014. [link](https://ieeexplore.ieee.org/document/6808537)


