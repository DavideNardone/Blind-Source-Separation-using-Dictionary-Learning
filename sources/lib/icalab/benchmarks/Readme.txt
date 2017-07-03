ICALAB for Signal Processing - BENCHMARKS            August 15, 2002

The directory Benchmarks contains typical signals and images 
for testing and comparison of various algorithms for 
ICA, BSE  and BSS.

The most interesting benchmarks are briefly described below.


1.  
ACsin10d.mat - contains 10 sine-wave sources of the following form:

s_n = sin ((2n-1)*w*k)  for  n = 1,2,....,10.

The sources can be easily separated by the second orders 
statistics methods (SOS) like AMUSE, EVD or SOBI algorithms. 
However, the higher order statistics (HOS)  
ICA algorithms fail to reconstruct such sources, because 
they are dependent.
It is interesting to note that different ICA algorithms 
(for example JADE and Natural gradient algorithms) give usually  
different (inconsistent) results ("independent components")
for this benchmark.

Acsin4d.mat benchmark is similar to  acsin10d.mat 
but it contains only 4 sources.


2.
ACsparse10.mat - contains 10 sparse (smooth bell-shape) sources 
that are approximately independent.  
The SOS blind source separation algorithms fail 
to separate such sources.
Also some ICA algorithms have failed to separate such sources.
Please try ICALAB to compare performance of various algorithms.


3.
ACvsparse10.mat - contains 10 very sparse (short regular pulses).
Second order BSS algorithms fail to separate these sources.


4. 
ABio7.mat - this benchmark contains 7 typical biological sources.
This benchmark was proposed by Allan Barros.


5. 
Sergio7.mat - this benchmark contains 7 sources (some of them 
are asymmetric distributed).  
This benchmark was proposed by Sergio Cruces.


6.
AC10-7sparse -contains 10 sensors signals which are mixture 
of 7 sources (extracted from the file ACsparse10.mat).


7. 
acspeech16.mat contains 16 typical speech signals which have 
temporal structure but they are not precisely independent.

Similar benchmarks are :

Speech4.mat, Speech8.mat. Speech10.mat and Speech20.mat with 
4, 8, 10 and 20 sounds (speech and music) sources.

10halo.mat contains 10 speech signals highly correlated 
(the all 10 speakers say the same sentence).


8. 
nband5.mat contains 5 narrow band sources.
This is standard "easy" benchamark.


9.
Gnband.mat contains 5 fourth order colored sources with distribution close 
to Gaussian. This is rather "difficult" benchmark. Please try program JADE-TD 
to separate the signals from their mixture.

   
10.
EEG19.mat consists 19 EEG signals with clear heart, eye movements 
and eye blinking artifacts.
