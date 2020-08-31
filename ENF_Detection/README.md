# Electric Network Frequency Detection in Audio Recordings

This package contains the MATLAB program for electric network frequency detection in audio recordings for audio forensic analysis (verified on MATLAB 2016b), including codes to reproduce figures showing the results of the [paper](https://ieeexplore.ieee.org/document/9143185).

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

## Citation Info:
  > \[1] G. Hua, H. Liao, Q. Wang, H. Zhang, and D. Ye, "Detection of electric network frequency in audio recordings – From theory to practical detectors," IEEE Trans. Inf. Forensics Security, DOI: 10.1109/TIFS.2020.3009579, 2020. [link](https://ieeexplore.ieee.org/document/9143185)
