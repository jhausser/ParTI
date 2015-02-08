%=====================================================================
% Programmer: 
% Tsung-Han Chan, E-mail: chantsunghan@gmail.com
% Date: May 20, 2009
% -------------------------------------------------------
% Reference: 
% T.-H. Chan, C.-Y. Chi, Y.-M. Huang and W.-K. Ma, ``A convex analysis 
% based minimum-volume enclosing simplex algorithm for hyperspectral 
% unmixing," to appear in IEEE Trans. Signal Processing, 2009. 
%======================================================================
% A implementation of MVES (Minimum-Volume Enclosing Simplex) Algorithm
% [A_est,S_est, time, iter_cnt] = MVES(X,N,show_flag)
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

function [A_est,S_est, time, iter_cnt] = MVES(X,N,show_flag)

t0 = clock;
%----------- Define default parameters------------------
TOL_obj = 1e-8; % convergence tolerance

%-------Step 1. and Step 2. Dimension reduction through affine set fitting-----------------
[M,L] = size(X);
d = mean(X,2);
U = X-d*ones(1,L);
R = U*U';
[eV D] = eig(R);
C = eV(:,M-N+2:end);
Xd = pinv(C)*(X-d*ones(1,L));  % Dimension reduced data vectors (Xd is (N-1)-by-L)

%---- Modification (finding boundary points using quickhull algorithm when N<10)--------
if N < 10
    Xi = convhulln(Xd');
    temp = vec(Xi);
    index = union(temp,temp(1));
    Lc = length(index);
    Xc = Xd(:,index);
else 
    Lc = L; Xc = Xd;
end

%--------Step 3. Find a feasible initialization to MVES-------------
c = [ones(Lc,1); sparse([],[],[],(N-1)*Lc,1,0)];
X_block = [Xc;-ones(1,Lc)];
A1 = kron(X_block',ones(1,N-1));
A2 = -kron(X_block',sparse(1:N-1,1:N-1,1,N-1,N-1));
AA = [A1;A2]';
K.l = N*Lc;
pars.fid = 0;
[x y] = sedumi(AA,0,c,K,pars);
H0 = reshape(y,N-1,N); % [H g]
H = H0(:,1:N-1);
g = H0(:,N);
obj0 = abs(det(H));

%----------Step 4. Solving the two linear programs--------------
i = 0; rec = 1; iter_cnt = 0;
while (rec > TOL_obj) & (iter_cnt<10*N)
    i = i+1;
    desire = rem(i,N-1)+(rem(i,N-1)==0)*(N-1);
    fixed = find((desire~=1:(N-1))>0);
    b = [];
    for j = 1:N-1;
        Hij = [H(1:desire-1,1:j-1),H(1:desire-1,j+1:N-1);H(desire+1:N-1,1:j-1),H(desire+1:N-1,j+1:N-1)];
        b = [b;(-1)^(desire+j)*det(Hij)];
    end
    b([b~=0]) = b([b~=0])./(max(abs(b([b~=0])))*ones(length(find(b~=0)),1)); 
    b = [b;0];
    H0_residual = H0(find((1:N-1~=desire)==1),:)';
    T = kron(ones(1,N-1),X_block');
    T1 = T(:,1:N);T2=T(:,N+1:end);
    c = [sparse([],[],[],Lc,1,0);ones(Lc,1)-T2*vec(H0_residual);0];
    A1 = [-X_block';T1;-b']';
    A2 = [-X_block';T1;b']';
    K.l = 2*Lc+1;
    [xa,ya] = sedumi(A1,b,c,K,pars);
    [xb,yb] = sedumi(A2,-b,c,K,pars);
    
    %------Step 5. ---------
    y = [ya yb];
    [val ind] = max([b'*ya -b'*yb]);
    y_true = y(:,ind);
    H0(desire,:) = y_true';
    H = H0(:,1:N-1); 
    
    %------Step 6.---------------------
    if desire == N-1
        iter_cnt = iter_cnt+1;
        rec = abs(obj0-abs(det(H)))/obj0;
        obj0 = abs(det(H));
        if show_flag
            disp(' ');
            disp(strcat('Number of iterations: ', num2str(iter_cnt)))
            disp(strcat('Relative change in objective: ', num2str(rec)))
            disp(strcat('Current volume: ', num2str(obj0)))
        end
    end
end

%--------Step. 7 to Step. 9------------------
g = H0(:,N);
alpha(:,N) = pinv(H)*g;
alpha(:,1:N-1) = pinv(H)+alpha(:,N)*ones(1,N-1);
A_est = C*alpha+d*ones(1,N);
S_est = [H*Xd-g*ones(1,L);ones(1,L)-ones(1,N-1)*(H*Xd-g*ones(1,L))];
%------------------------------------------
time = etime(clock,t0);