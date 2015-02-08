MVES (Minimum-volume Enclosing Simplex) algorithm ReadMe

Our MVES algorithm relies on a (also free) convex optimization software SeDuMi. 
Please install it on your computer if you don't have one.

The main function for the MVES algorithm is 

[A_est,S_est, time, iter_cnt] = MVES(X,N,show_flag)

%======================================================================
%  Input
%  X is M-by-L data matrix where M is the spectral bands (or observations) and L is the number of pixels (data length).   
%  N is the number of endmembers (or sources).
%  show_flag: 1- display current information in MVES, and 0 - otherwise 
%----------------------------------------------------------------------
%  Output
%  A_est is M-by-N: estimated endmember signatures (or mixing matrix) obtained by MVES.
%  S_est is N-by-L: estimated abundance vectors (or source matrix) obtained by MVES.
%  time is the computation time (in secs). 
%  iter_cnt is the passed number of iterations in MVES. 
%========================================================================

Any questions and suggestions are very welcome. Please send to 
Tsung-Han Chan, email: chanchantsunghan@gmail.com

