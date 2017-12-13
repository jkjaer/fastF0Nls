# fastF0Nls - fast nonlinear least squares estimation of the fundamental frequency (aka pitch)
Periodic signals are encountered in many real-world applications such as music processing, speech processing, sonar, order analysis, and electrocardiography (ECG). Such signals can be modelled as a weighted sum of sinusoids whose frequencies are integer multiples of a common fundamental frequency which in audio and speech applications is often referred to as the pitch. Therefore, an important and fundamental problem in the above mentioned applications is to estimate this fundamental frequency from an observed (and often noisy) data set.

Multiple estimation methods have been proposed in the scientific literature ranging from simple correlation-based methods to parametric methods. Although the parametric methods in general are much more accurate and robust to noise than the correlation-based methods, they suffer from a high computational complexity. Consequently, the correlation-based methods remain very popular despite that they require all sorts of heuristic post-processing to give a satisfactory performance.

The nonlinear least squares (NLS) estimator of the fundamental frequency has been known since 1992 and been shown to one of the most accurate and robust fundamental frequency estimator, in particular for short segments of data. This git repository contains C++- and MATLAB-implementations of a very fast algorithm for computing the NLS estimate and the number of harmonic component. The algorithm is described in 
[1]

[2]
