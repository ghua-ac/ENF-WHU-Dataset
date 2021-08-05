# About

This repo contains the ENF-WHU audio recording dataset collected around Wuhan University campus and the MATLAB programs for electronic network frequency (ENF) detection, enhancement, and robust estimation, in ENF-based audio forensic applications.

## ENF-WHU Dataset
- **Recording location:** classroom, campus path, meeting room, graduate student office, dormitory, library.
- **Environment diversity:** day/night, rainy/suny, interior/exterior.
- **Recording device:** popular smartphone and voice recorder.
- **Duration:** 5~20 minutes
- **Format:** PCM WAVE
- **Quantization depth:** 16-bit
- **Channel:** mono
- **Sampling frequencuy:** 8000 Hz (400 Hz for reference data)
- **Category:**<br>
  **H1:** "001\~130.wav" 130 real-world recordings with captured (noisy) ENF.<br>
  **H1_ref:** "001_ref\~130_ref.wav" the corresponding 130 reference ENF (noise-free) obtained from power main.<br>
  **H0:** "O01\~O10.wav" 10 real-world recordings without captured ENF. "01\~40.wav" 40 segments under H0 obtained by random cropping the 10 recordings.
  **H1_ref_one_day:** the corresponding one-day (24 hours) reference ENF for the 130 recordings. "003-004_ref.wav" means "003.wav" and "004.wav" in **H1** are recorded within the same day.

## MATLAB Programs
### ENF Detection
- Clairvoyant detectors: NP detectors assuming perfect knowledge of ENF.
	1. GMF: a standard NP detector.
	2. MF-like approximation: avoid the requirement of unknown noise covariance matrix.
	3. Asymptotic approximation: trade-off between computational complexity and detection performance.
- GLRT detectors: ENF assumed unknown and deterministic.
	1. LS-LRT: MF-like with unknown parameters replaced by the MLEs.
	2. naive-LRT: MF-like with the unknown IFs replaced by nominal value.
- TF domain detector: ENF assumed unknown and random.
	- Test statistic is the sample variance of the strongest time-frequency line (e.g., STFT + peak)
	- Exploiting slow-varying nature of ENF, thus test statistic is large under H0 and small under H1. 

### ENF Enhancement and Estimation

It contains our proposed ENF enhancement and estimation methods including
- Single-tone model based ENF enhancement method [3],
- Multi-tone harmonic model based enhancement and harmonic selection for robust ENF estimation [2],<br>
in comparison with the following existing works 
- Adaptive multi-trace carving (AMTC) method by [Q. Zhu <em>et al.</em> 2021 IEEE TIFS](https://ieeexplore.ieee.org/document/9220114),
- Multi-tone maximum likelihood estimator (MLE) by [Bykhovsky and Cohen 2013 IEEE TIFS](https://ieeexplore.ieee.org/document/6482617),
- Multi-tone weighted maximum likelihood estimator (WMLE) by [Hajj-Ahmad <em>et al.</em> 2013 IEEE SPL](https://ieeexplore.ieee.org/document/6557080),<br>
evaluated using both synthetic data and the real-world recordings from the ENF-WHU dataset.

## Citation Information
- **ENF Detection:**
 > \[1] G. Hua, H. Liao, Q. Wang, H. Zhang, and D. Ye, "Detection of electric network frequency in audio recordings – From theory to practical detectors," IEEE Trans. Inf. Forensics Security, vol. 16, pp. 236–248, 2021. [link](https://ieeexplore.ieee.org/document/9143185)
- **ENF Enhancement:**
 > \[2] G. Hua, H. Liao, H. Zhang, D. Ye, and J. Ma, "Robust ENF estimation based on harmonic enhancement and maximum weight clique," IEEE Trans. Inf. Forensics Security, DOI: 10.1109/TIFS.2021.3099697, 2021. [link](https://ieeexplore.ieee.org/document/9494518)<br>
 > \[3] G. Hua and H. Zhang, "ENF signal enhancement in audio recordings," IEEE Trans. Inf. Forensics Security, vol. 15, pp. 1868-1878, 2020. [link](https://ieeexplore.ieee.org/document/8894138)
- **Related Works:**
 > \[4] G. Hua, "Error analysis of forensic ENF matching," in Proc. 2018 IEEE International Workshop on Information Forensics and Security (WIFS), pp. 1-7, Hong Kong, Dec. 2018. [link](https://ieeexplore.ieee.org/document/8630786)<br>
 > \[5] G. Hua, G. Bi, and V. L. L. Thing, "On practical issues of electric network frequency based audio forensics," IEEE Access, vol. 5, pp. 20640-20651, Oct. 2017. [link](https://ieeexplore.ieee.org/document/7807225)<br>
 > \[6] G. Hua, Y. Zhang, J. Goh, and V. L. L. Thing, "Audio authentication by exploring the absolute error map of the ENF signals," IEEE Trans. Inf. Forensics Security, vol. 11, no. 5, pp. 1003-1016, May 2016. [link](https://ieeexplore.ieee.org/document/7378470)<br>
 > \[7] G. Hua, J. Goh, and V. L. L. Thing, “A dynamic matching algorithm for audio timestamp identification using the ENF criterion,” IEEE Trans. Inf. Forensics Security, vol. 9, no. 7, pp. 1045-1055, Jul. 2014. [link](https://ieeexplore.ieee.org/document/6808537)


