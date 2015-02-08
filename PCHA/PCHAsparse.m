function [XC,S,C,SSE,varexpl]=PCHAsparse(X,noc,varargin)
% Principal Convex Hull Analysis (PCHA) / Archetypal Analysis for sparse
% data where zero entries are treated as missing (useful for matrix completion such as the Netflix and Movielens problem)
%
% Written by Morten Mørup
%
% Usage:
%   [C,S,SSE,varexpl]=PCHAsparse(X,noc,varargin)
%
%   Solves the following PCHA/AA problem
%   \sum_(X)_ij>0 [(X)_ij-(XCS)_ij]^2 s.t. |s_j|_1=1, |c_j|_1=1, S>=0, C>=0
%   treating entries in X that are zero as missing
%
% Input:
% X             data array (Missing entries set to zero)
% noc           number of components
% opts.         Struct containing:
%       C            initial solution (optional) (see also output)
%       S            initial solution (optional) (see also output)
%       maxiter      maximum number of iterations (default: 500 iterations)
%       conv_crit    The convergence criteria (default: 10^-6 relative change in costfcn)
%
% Output:
% XC            I x noc feature matrix (i.e. XC=X*C forming the archetypes) 
% S             noc x J matrix, S>=0 |S_j|_1=1
% C             J x noc matrix, C>=0 |C_j|_1=1
% SSE           Sum of Squares Error
% varexpl       Percent variation explained by the model
%
% Copyright (C) Morten Mørup and Technical University of Denmark, 2010
warning('off','MATLAB:dispatcher:InexactMatch')

[I,J]=size(X);
nmissing=I*J-nnz(X);
SST=sum(sum(X.*X));
[Ix,Jx,valx]=find(X);
M=sparse(Ix,Jx,ones(size(valx)),I,J);

% Set algorithm parameters
if nargin>=3, opts = varargin{1}; else opts = struct; end
conv_crit=mgetopt(opts,'conv_crit',10^-6);
maxiter=mgetopt(opts,'maxiter',500);

% Initilize C and S
if isfield(opts,'C')    
    C = opts.C;
    MC=M*C+1e-3;
    XC=(X*C)./MC;
else      
   C=-log(rand(J,noc));
   C=C./repmat(sum(C),J,1);
   MC=M*C+1e-3;
   XC=(X*C)./MC;
end
if isfield(opts,'S') 
    S=full(opts.S);
else   
    S=-log(rand(noc,J));
    S=S./repmat(sum(S),noc,1);
end

% Set PCHA parameters
iter=0;
dSSE=inf;
Rec=dprod(M,XC,S);                       
SSE=full(SST+sum(sum((Rec-2*X).*Rec)));   
t1=cputime;
muC=1/(SST/J);
muS=muC;
varexpl=0;

% Display algorithm profile
disp([' '])
disp(['Principal Conical Hull analysis'])
disp(['A ' num2str(noc) ' component model will be fitted']);
disp([ num2str(nmissing/(I*J)*100) ' pct. of data is missing '])
disp(['To stop algorithm press control C'])
disp([' ']);
dheader = sprintf('%12s | %12s | %12s | %12s | %12s | %12s | %12s','Iteration','Expl. var.','Cost func.','Delta costf.','muC','muS',' Time(s)   ');
dline = sprintf('-------------+--------------+--------------+--------------+--------------+--------------+--------------+');


while abs(dSSE)>=conv_crit*abs(SSE) && iter<maxiter && varexpl<0.9999
    if mod(iter,10)==0
         disp(dline); disp(dheader); disp(dline);
    end
    told=t1;
    iter=iter+1;
    SSE_old=SSE;
        
    % C update    
    [C,SSE,muC,Rec,XC,MC]=Cupdate(M,C,X,S,XC,MC,Rec,muC,SST,SSE,1);
    
    % S update    
    [S,SSE,muS,Rec]=Supdate(M,S,X,Rec,XC,muS,SST,SSE,1);
    
    
    % Evaluate and display iteration
    dSSE=SSE_old-SSE;
    t1=cputime;
    if rem(iter,1)==0
        varexpl=(SST-SSE)/SST;
        fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4e | %12.4f \n',iter,varexpl,SSE,dSSE/abs(SSE),muC,muS,t1-told);
    end
end

% display final iteration
varexpl=(SST-SSE)/SST;
disp(dline);
disp(dline);
fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4e | %12.4f \n',iter,varexpl,SSE,dSSE/abs(SSE),muC,muS,t1-told);

% sort components according to importance
[val,ind]=sort(sum(S,2),'descend');
S=S(ind,:);
C=C(:,ind);

% -------------------------------------------------------------------------
% Parser for optional arguments
function var = mgetopt(opts, varname, default, varargin)
if isfield(opts, varname)
    var = getfield(opts, varname); 
else
    var = default;
end
for narg = 1:2:nargin-4
    cmd = varargin{narg};
    arg = varargin{narg+1};
    switch cmd
        case 'instrset',
            if ~any(strcmp(arg, var))
                fprintf(['Wrong argument %s = ''%s'' - ', ...
                    'Using default : %s = ''%s''\n'], ...
                    varname, var, varname, default);
                var = default;
            end
        otherwise,
            error('Wrong option: %s.', cmd);
    end
end


% -------------------------------------------------------------------------
function [S,SSE,muS,Rec]=Supdate(M,S,X,Rec,XC,muS,SST,SSE,niter)
    noc=size(S,1);
    for k=1:niter
        SSE_old=SSE;
        g=XC'*(Rec-X);       
        g=g-repmat(sum(g.*S),noc,1);
        stop=0;
        Sold=S;
        while ~stop
            S=Sold-g*muS;
            S(S<0)=0;
            S=S./repmat(sum(S),noc,1);
            Rec=dprod(M,XC,S);                       
            SSE=full(SST+sum(sum((Rec-2*X).*Rec)));            
            if SSE<=SSE_old*(1+1e-9)
                muS=muS*1.2;
                stop=1;
            else
                muS=muS/2;
            end
        end
    end
    
%--------------------------------------------------------------------
function [C,SSE,muC,Rec,XC,MC]=Cupdate(M,C,X,S,XC,MC,Rec,muC,SST,SSE,niter)
    J=size(C,1);
    for k=1:niter
        SSE_old=SSE;
        T=(Rec-X)*S';
        g=X'*(T./MC)-M'*(T.*XC./MC);
        g=g-repmat(sum(g.*C),J,1);
        stop=0;
        Cold=C;
        while ~stop            
            C=Cold-muC*g;
            C(C<0)=0;
            C=C./repmat(sum(C),J,1);            
            MC=(M*C+1e-3);
            XC=(X*C)./MC;            
            Rec=dprod(M,XC,S);
            SSE=full(SST+sum(sum((Rec-2*X).*Rec)));            
            if SSE<=SSE_old*(1+1e-9)
                muC=muC*1.2;
                stop=1;
            else
                muC=muC/2;                
            end
        end
    end
    
%--------------------------------------------------------------------
function R=dprod(Q,W,H)
% Direct product
W=W';
[I,J,val]=find(Q);
step=10000;
for k=1:ceil((length(I)/step))
    ind=(k-1)*step+1:min([k*step, length(I)]);
    val(ind)=sum(W(:,I(ind)).*H(:,J(ind)));
end
R=sparse(I,J,val,size(Q,1),size(Q,2)); 




