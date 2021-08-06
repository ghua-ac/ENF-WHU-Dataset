### ENF Enhancement and Estimation

It contains our proposed ENF enhancement and estimation methods including
- Single-tone model based ENF enhancement method [3],
- Multi-tone harmonic model based enhancement and harmonic selection for robust ENF estimation [2],<br>

in comparison with the following existing works 

- Adaptive multi-trace carving (AMTC) method by [Q. Zhu <em>et al.</em> 2021 IEEE TIFS](https://ieeexplore.ieee.org/document/9220114),
- Multi-tone maximum likelihood estimator (MLE) by [Bykhovsky and Cohen 2013 IEEE TIFS](https://ieeexplore.ieee.org/document/6482617),
- Multi-tone weighted maximum likelihood estimator (WMLE) by [Hajj-Ahmad <em>et al.</em> 2013 IEEE SPL](https://ieeexplore.ieee.org/document/6557080),<br>

evaluated using both synthetic data and the real-world recordings from the ENF-WHU dataset.

### Usage

- The MATLAB implementation of Bron-Kerbosch algorithm for maximal clique detection
```
func_BK_MaxClique.m
func_BK_MaxIS.m
```
is incorporated from the work of [Berk Birand](https://kr.mathworks.com/matlabcentral/fileexchange/24591-bron-kerbosch-maximal-independent-set-and-maximal-clique-algorithms)


### Citation Info:
 > \[1] G. Hua, H. Liao, H. Zhang, D. Ye, and J. Ma, "Robust ENF estimation based on harmonic enhancement and maximum weight clique," IEEE Trans. Inf. Forensics Security, DOI: 10.1109/TIFS.2021.3099697, 2021. [link](https://ieeexplore.ieee.org/document/9494518)<br>
 > \[2] G. Hua and H. Zhang, "ENF signal enhancement in audio recordings," IEEE Trans. Inf. Forensics Security, vol. 15, pp. 1868-1878, 2020. [link](https://ieeexplore.ieee.org/document/8894138)
