function [A_est time]=SDVMM(Y,N,r)     
%=====================================================================
% Programmers: 
% Tsung-Han Chan, E-mail: thchan@ieee.org
% Ji-Yuan Liou
% Date: Aug., 2011
%======================================================================
% A implementation of WCR-SVMAX
% [A_est time]=SDVMM(Y,N,r)
%======================================================================
%  Input
%  Y is M-by-L data matrix where M is the spectral bands (or observations) and L is the number of pixels (data length).   
%  N is the number of endmembers (or sources).
%  r is a back-off tolerance which is given by r = 1.3*(noise standard deviation).
%----------------------------------------------------------------------
%  Output
%  A_est is M-by-N: estimated endmember signatures (or mixing matrix) obtained by SDVMM.
%  time is the computation time (in secs). 
%========================================================================

t0 = clock;
[M,L] = size(Y);
d = mean(Y,2);
U = Y-d*ones(1,L);
OPTS.disp = 0;
[C D] = eigs(U*U',N-1,'LM',OPTS);
Yd_t = C'*U;                                                % dimension reduced data
% ========= SDVMM algorithm =========
H_set=[];  index = []; P = eye(N); Yd = [Yd_t; ones(1,L)];  % step S1 
for i=1:N   
    norm_vectors = sum(abs(P*Yd).^2).^(1/2);                % step S2 or S4
    fesi_index = find( norm_vectors > r );                  % fesi_index represents the feasible index set Nj in step S2 or S4
    [val ind_p] = max(norm_vectors(:,fesi_index));
    ind = fesi_index(ind_p); 
    z = (r/val)*(P*Yd(:,ind));              
    z(N) = 0;                                               % step S3 or S5
    H_set = [H_set Yd(:,ind)-z];                             
    P = eye(N) - H_set*pinv(H_set);                                
    index = [index ind];                                        
end
A_est=C*H_set(1:N-1,:)+d*ones(1,N);
time = etime(clock,t0);

