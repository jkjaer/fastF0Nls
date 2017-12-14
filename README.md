# fastF0Nls - fast nonlinear least squares estimation of the fundamental frequency (aka pitch)

Periodic signals are encountered in many real-world applications such as music processing, speech processing, sonar, order analysis, and electrocardiography (ECG). Such signals can be modelled as a weighted sum of sinusoids whose frequencies are integer multiples of a common fundamental frequency which in audio and speech applications is often referred to as the pitch. Therefore, an important and fundamental problem in the above mentioned applications is to estimate this fundamental frequency from an observed (and often noisy) data set.

Multiple estimation methods have been proposed in the scientific literature ranging from simple correlation-based methods (such as Yin, PEFAC, and RAPT) to parametric methods (such as subspace methods, filtering methods, harmonic summation, and NLS). Although the parametric methods in general are much more accurate and robust to noise than the correlation-based methods, they typically suffer from a high computational complexity. Consequently, the correlation-based methods remain very popular despite that they require all sorts of heuristic post-processing to give a satisfactory performance.

The nonlinear least squares (NLS) estimator of the fundamental frequency is a parametric method which has been known since 1992 and been shown to be one of the most accurate and robust fundamental frequency estimators, in particular for short segments of data. Unfortunately, the computational complexity of any published NLS algorithm is very high so most people have instead used harmonic summation which is an approximate NLS method. Recently, however, we published a very fast algorithm for computing the NLS estimate and the number of harmonic components, and this git repository contains C++- and MATLAB-implementations of this algorithm. The details of the algorithm can be found in

- Fast fundamental frequency estimation: Making a statistically efficient estimator computationally efficient. Nielsen, Jesper Kjær; Jensen, Tobias Lindstrøm; Jensen, Jesper Rindom; Christensen, Mads Græsbøll; Jensen, Søren Holdt. In: Signal Processing, 135, 2017, pp. 188-197.
- Bayesian Model Comparison With the g-Prior. Nielsen, Jesper Kjær; Christensen, Mads Græsbøll; Cemgil, Ali Taylan; Jensen, Søren Holdt. In: IEEE Transactions on Signal Processing, 62 (1), 2014, pp. 225-238.

In addition to these papers, we have also made a one hour video tutorial on fundamental frequency in general and our fast algorithm in particular. The video is published on [YouTube](https://www.youtube.com/watch?v=F0XgU-9ERp4).

To get started with the code, please see the examples and the documentation.

## Scope and limitations

Please note that the code only contains a pitch estimator and NOT a pitch tracker. The difference is that a tracker contains a smoothing step on top of the estimator. The smooting step is there to minimise the risk of, e.g., octave errors (aka pitch halving) by smoothing out the estimates produced by the estimator which typically analyse the data on a segment-by-segment basis. Of course, our estimator can also be used inside a pitch tracker, and we recommend that the smoothing step by Tabrikian et al. is used. See more in

- Maximum a-posteriori probability pitch tracking in noisy environments using harmonic model, Tabrikian, Joseph; Dubnov, Shlomo; and Dickalov, Yulya. In: IEEE Transactions on Speech and Audio Processing 12 (1), 2004, pp. 76-87.

For white Gaussian noise (WGN), the NLS estimator is the maximum likelihood estimator and is, therefore, asymptotically optimal (in a statistical sense). That is, no other unbiased estimator of the fundamental frequency has a lower variance than the NLS estimator if enough data are observed and the noise is white and Gaussian. In our experience, the NLS estimator is also one of the best estimators for short data segments, and it works well even for situations where a data segment contain only one cycle of the periodic signal. For voiced speech, where the lowest fundamental frequency is typically bigger than 80 Hz, this means that the estimator typically works well down to a segment length of 12.5 ms. Increasing the segment length will of course increase the estimation accuracy and the robustness to noise, provided that the signal is approximately stationary.

In our experience, the estimator does typically not break down if the noise is not white and Gaussian. However, if the noise is coloured and has most of the energy at the lower frequencies, then the estimator can suffer from problems with octave errors. In this case, we recommend that some kind of pre-whitening is applied to the data prior to estimating the fundamental frequency.

## License

The code is published under a BSD 3-clause "New" or "Revised" License. This is a very permissive license that enables the code to be used free of charge in both academia and industry. Besides honouring the license, we therefore also expect that you 

- cite our work if you use the code academically, and
- notify us if you use the code commercially.

Constructive feedback in any form is also very much appreciated. This applies to everything from bugs to data where the estimator fails. Please send your feedback to Jesper (jkn@create.aau.dk).

The code has been written and is maintained by

- Tobias Lindstrøm Jensen, Intel Mobile Communications Denmark, Aalborg, Denmark
- Jesper Kjær Nielsen, Audio Analysis Lab, CREATE, Aalborg University, Denmark
