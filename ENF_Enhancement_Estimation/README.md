### ENF Enhancement and Estimation

It contains our proposed ENF enhancement and estimation methods including
- Single-tone model based ENF enhancement method [3],
- Multi-tone harmonic model based enhancement and harmonic selection for robust ENF estimation [2],<br>

in comparison with the following existing works 

- Adaptive multi-trace carving (AMTC) method by [Q. Zhu <em>et al.</em> 2021 IEEE TIFS](https://ieeexplore.ieee.org/document/9220114),
- Multi-tone maximum likelihood estimator (MLE) by [Bykhovsky and Cohen 2013 IEEE TIFS](https://ieeexplore.ieee.org/document/6482617),
- Multi-tone weighted maximum likelihood estimator (WMLE) by [Hajj-Ahmad <em>et al.</em> 2013 IEEE SPL](https://ieeexplore.ieee.org/document/6557080),

evaluated using both synthetic data and the real-world recordings from the ENF-WHU dataset.

### Usage

The MATLAB implementation of Bron-Kerbosch algorithm for maximal clique detection
 > func_BK_MaxClique.m<br>
 > func_BK_MaxIS.m<br>

is incorporated from the work of [Berk Birand](https://www.mathworks.com/matlabcentral/fileexchange/24591-bron-kerbosch-maximal-independent-set-and-maximal-clique-algorithms).
 > func_BPF.p

creates the comb filter for harmonic bandpass filtering centered at 100, 150, 200, 250, 300, 350 Hz.
 > func_ENF_synthesis_corrupted_harmonic.m

creates a synthetic time-domain ENF waveform with a few corrupted harmonics.
 > func_RFA.p
 > func_RFA_multi.p

are the implementations of the proposed robust filtering algorithm (RFA) [2] and harmonic robust filtering algorithm (HRFA) [1].
 > func_STFT_multi_tone_MWC.p

combines the HRFA and graph-based harmonic selection algorithm (GHSA), which gives the full implementation of the framework in [1].
 > func_STFT_multi_tone_component.p
 
calculates the harmonic components respectively.
 > func_STFT_multi_tone_search.m<br>
 > func_STFT_multi_tone_search_weighted.m

are implementations the MLE [Bykhovsky and Cohen 2013 IEEE TIFS](https://ieeexplore.ieee.org/document/6482617) and WMLE [Hajj-Ahmad <em>et al.</em> 2013 IEEE SPL](https://ieeexplore.ieee.org/document/6557080), respectively.
 > func_STFT_single_tone.m
 
 is the implementation of single-tone MLE.
 > func_empirical_coeff_thre.m

calculates the empirical threshold &eta;<sub>R</sub> in (20) of [1].

### Citation Info:
 > \[1] G. Hua, H. Liao, H. Zhang, D. Ye, and J. Ma, "Robust ENF estimation based on harmonic enhancement and maximum weight clique," IEEE Trans. Inf. Forensics Security, DOI: 10.1109/TIFS.2021.3099697, 2021. [link](https://ieeexplore.ieee.org/document/9494518)<br>
 > \[2] G. Hua and H. Zhang, "ENF signal enhancement in audio recordings," IEEE Trans. Inf. Forensics Security, vol. 15, pp. 1868-1878, 2020. [link](https://ieeexplore.ieee.org/document/8894138)
