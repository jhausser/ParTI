function [A_est time]=ADVMM(Y,N,r,show_flag)
%=====================================================================
% Programmers: 
% Tsung-Han Chan, E-mail: thchan@ieee.org  
% Ji-Yuan Liou
% Date: Aug., 2011
%======================================================================
% A implementation of ADVMM
% [A_est,time] = ADVMM(Y,N,r,show_flag)
%======================================================================
%  Input
%  Y is M-by-L data matrix where M is the spectral bands (or observations) and L is the number of pixels (data length).   
%  N is the number of endmembers (or sources).
%  r is a back-off tolerance which is given by r = 1.3*(noise standard deviation).
%  show_flag: 1- display current information in ADVMM, and 0 - otherwise 
%----------------------------------------------------------------------
%  Output
%  A_est is M-by-N: estimated endmember signatures (or mixing matrix) obtained by ADVMM.
%  time is the computation time (in secs). 
%========================================================================

TOL_CONVER =5e-5;                               % convergence tolerence
t0 = clock;
[M,L] = size(Y);
d = mean(Y,2);
U_obs = Y - d*ones(1,L);
OPTS.disp = 0;
[C D] = eigs(U_obs*U_obs',N-1,'LM',OPTS);
Y_tilde = C'*U_obs;                             % dimension reduced data 
% =========================== ADVMM algorithm =========================
%-----------------------------step S1: initialization-------------------
order = randperm(L);
B = order(1,1:N);
V = Y_tilde(:,B); 
U = zeros(N-1,N);       
H = [V-U;ones(1,N)];
%-----------------------------step S2-----------------------------------
volume = abs(det(H));

y = 0;   rec = 1;  iter_cnt = 0;   
while (rec > TOL_CONVER) & (iter_cnt < 2*N)    % convergence criterion (step S5)
    y = y + 1;
    j = rem(y,N) + (rem(y,N)==0)*N;   % j = 1 ~ N
    kj = [];
    %--------------------step S3----------------------
    for i = 1:N-1
        Qij = [H(1:i-1,1:j-1),H(1:i-1,j+1:N);H(i+1:N,1:j-1),H(i+1:N,j+1:N)];
        kj = [kj;(-1)^(j+i)*det(Qij)];          % compute kj
    end
    U(:,j) = (r/norm(kj))*kj;           
    [maxvalue max_index] = max(kj'*Y_tilde);
    V(:,j) = Y_tilde(:,max_index);    
    H=[V-U;ones(1,N)];                          % update
    %--------------------step S4----------------------
    if j == N
        iter_cnt = iter_cnt+1;                  % the no. of outer iteration 
        rec = abs(volume-abs(det(H)))/volume;   % relative change in determinant value
        volume = abs(det(H));      
        if show_flag
            disp(' ');
            disp(strcat('Number of iterations: ', num2str(iter_cnt)));
            disp(strcat('Relative change in objective: ', num2str(rec)));
            disp(strcat('Current volume: ', num2str(volume)));
        end
    end
end
alpha_est = V - U;
A_est = C*alpha_est + d*ones(1,N);
time = etime(clock,t0);
