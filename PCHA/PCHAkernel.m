function [S,C,SSE,varexpl]=PCHAkernel(K,noc,I,U,delta,varargin)
% Kernel Principal Convex Hull Analysis (PCHA) / Kernel Archetypal Analysis
%
% Written by Morten Mørup
%
% Usage:
%   [S,C,SSE,varexpl]=PCHAkernel(K,noc,I,U,delta,opts)
%
% Input:
% K             Kernel Matrix
% noc           number of components
% I             Entries of K to use for dictionary in C (default: I=1:size(K,2))
% U             Entries of K to model in S              (default: U=1:size(K,2))
% delta         relaxation of C, i.e. 1-delta<=|C_j|_1<=1+delta (default: delta=0, i.e. |C_j|_1=1)
% opts.         Struct containing:
%       C            initial solution (optional) (see also output)
%       S            initial solution (optional) (see also output)
%       maxiter      maximum number of iterations (default: 500 iterations)
%       conv_crit    The convergence criteria (default: 10^-6 relative change in costfcn)
%
% Output:
% S             noc x length(U) matrix, S>=0 |S_j|_1=1
% C             length(I) x noc matrix, C>=0 1-delta<=|C_j|_1<=1+delta
% SSE           Sum of Squares Error
% varexpl       Percent variation explained by the model
%
% Copyright (C) Morten Mørup and Technical University of Denmark, 2010

warning('off','MATLAB:dispatcher:InexactMatch')
% Remove missing entries from X

if nargin<5
    delta=0;
end
if nargin<4
    U=1:size(K,2);
end
if nargin<3
    I=1:size(K,2);
end

% Set algorithm parameters
if nargin<5
    delta=0;
end
if nargin>=6, opts = varargin{1}; else opts = struct; end
conv_crit=mgetopt(opts,'conv_crit',10^-6);
maxiter=mgetopt(opts,'maxiter',500);
SST=trace(K(U,U));

% Initilize C and S
if isfield(opts,'C')
    C = opts.C;
else
    % Initialize by FurthestSum
    exclude=setdiff(1:size(K,2),I);
    i=FurthestSum(K,noc,ceil(size(K,2)*rand),exclude);
    C=zeros(length(I),noc);
    for k=1:noc
        ind=find(i(k)==I);
        C(ind,k)=1;
    end        
end
KC=K(:,I)*C; 
muS=1;
muC=1;
mualpha=1;

if isfield(opts,'S')
    S=opts.S;      
    CtKC=C'*K(I,I)*C;    
    SSt=S*S';           
    SSE=SST-2*sum(S'.*KC(U,:))+sum(sum(CtKC.*SSt));                
else       
    CtKC=C'*K(I,I)*C;    
    S=-log(rand(noc,length(U)));
    S=S./repmat(sum(S),noc,1);
    SSt=S*S';
    SSE=SST-2*sum(sum(S'.*KC(U,:)))+sum(sum(CtKC.*SSt));                
    [S,SSE,muS,SSt]=Supdate(S,KC(U,:),CtKC,muS,SST,SSE,25);         
end

% Set PCHA parameters
iter=0;
dSSE=inf;
t1=cputime;
varexpl=(SST-SSE)/SST;

% Display algorithm profile
disp([' '])
disp(['Principal Conical Hull analysis'])
disp(['A ' num2str(noc) ' component model will be fitted']);
disp(['To stop algorithm press control C'])
disp([' ']);
dheader = sprintf('%12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s','Iteration','Expl. var.','Cost func.','Delta costf.','muC','muAlpha','muS',' Time(s)   ');
dline = sprintf('-------------+--------------+--------------+--------------+--------------+--------------+--------------+--------------+');


while abs(dSSE)>=conv_crit*abs(SSE) && iter<maxiter  && varexpl<0.9999
    if mod(iter,10)==0
         disp(dline); disp(dheader); disp(dline);
    end
    told=t1;
    iter=iter+1;
    SSE_old=SSE;
        
    % C update    
    [C,SSE,muC,mualpha,CtKC,KC]=Cupdate(K,I,U,S,KC,SSt,C,delta,muC,mualpha,SST,SSE,5);
    
    % S update        
    [S,SSE,muS,SSt]=Supdate(S,KC(U,:),CtKC,muS,SST,SSE,5);            
    
    % Evaluate and display iteration
    dSSE=SSE_old-SSE;
    t1=cputime;
    if rem(iter,1)==0
        varexpl=(SST-SSE)/SST;
        fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4e | %12.4e | %12.4f \n',iter,varexpl,SSE,dSSE/abs(SSE),muC,mualpha,muS,t1-told);
    end
end

% display final iteration
varexpl=(SST-SSE)/SST;
disp(dline);
disp(dline);
fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4e | %12.4e | %12.4f \n',iter,varexpl,SSE,dSSE/abs(SSE),muC,mualpha,muS,t1-told);

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
function [S,SSE,muS,SSt]=Supdate(S,KC,CtKC,muS,SST,SSE,niter)
    [noc,J]=size(S);
    e=ones(noc,1);    
    for k=1:niter
        SSE_old=SSE;
        g=(CtKC*S-KC')/(SST/J);         
        g=g-e*sum(g.*S);        
        stop=0;
        Sold=S;
        while ~stop
            S=Sold-g*muS;
            S(S<0)=0;
            S=S./repmat(sum(S),noc,1);            
            SSt=S*S';          
            SSE=SST-2*sum(sum(S'.*KC))+sum(sum(CtKC.*SSt));            
            if SSE<=SSE_old*(1+1e-9)
                muS=muS*1.2;
                stop=1;
            else
                muS=muS/2;
            end
        end
    end

%--------------------------------------------------------------------
function [C,SSE,muC,mualpha,CtKC,KC]=Cupdate(K,I,U,S,KC,SSt,C,delta,muC,mualpha,SST,SSE,niter)
    JJ=size(K,2);
    [J,noc]=size(C);
    if nargin<12
        niter=1;
    end          
    if delta~=0
        alphaC=sum(C);    
        C=C*diag(1./alphaC);
    end
    e=ones(J,1);
    if length(U)~=JJ
        KSt=K(:,U)*S';  
    else      
        KSt=K*S';
    end
    for k=1:niter 
        
        % Update C        
        SSE_old=SSE;        
        if length(I)~=JJ
            g=(KC(I,:)*SSt-KSt(I,:))/SST;       
        else
            g=(KC*SSt-KSt)/SST;       
        end
        if delta~=0
            g=g*diag(alphaC);
        end
        g=g-e*sum(g.*C);        
        stop=0;
        Cold=C;
        while ~stop
            C=Cold-muC*g;
            C(C<0)=0;
            
            nC=sum(C)+eps;            
            C=C*sparse(1:noc,1:noc,1./(nC));
            if delta~=0
                Ct=sparse(C*diag(alphaC));
            else
                Ct=sparse(C);
            end
            if length(I)~=JJ
                KC=K(:,I)*Ct;
                CtKC=Ct'*KC(I,:);
            else
                KC=K*Ct;
                CtKC=Ct'*KC;
            end
            if length(U)~=JJ
                SSE=SST-2*sum(sum(S'.*KC(U,:)))+sum(sum(CtKC.*SSt));                                    
            else
                SSE=SST-2*sum(sum(S'.*KC))+sum(sum(CtKC.*SSt));                                    
            end
            if SSE<=SSE_old*(1+1e-9)
                muC=muC*1.2;
                stop=1;
            else
                muC=muC/2;
            end
        end        
        
        % Update alphaC        
        SSE_old=SSE;
        if delta~=0   
            if length(I)~=JJ
                g=(sum(CtKC.*SSt)./alphaC-sum(C.*KSt(I,:)))/(SST*J);                       
            else
                g=(sum(CtKC.*SSt)./alphaC-sum(C.*KSt))/(SST*J);                       
            end
            stop=0;
            alphaCold=alphaC;
            while ~stop
                alphaC=alphaCold-mualpha*g;
                alphaC(alphaC<1-delta)=1-delta;
                alphaC(alphaC>1+delta)=1+delta;                            
                KCt=KC*diag(alphaC./alphaCold);
                CtKCt=diag(alphaC./alphaCold)*CtKC*diag(alphaC./alphaCold);            
                if length(U)~=JJ
                    SSE=SST-2*sum(sum(S'.*KCt(U,:)))+sum(sum(CtKCt.*SSt));                
                else
                    SSE=SST-2*sum(sum(S'.*KCt))+sum(sum(CtKCt.*SSt));                
                end
                if SSE<=SSE_old*(1+1e-9)  
                    mualpha=mualpha*1.2;
                    stop=1;
                    KC=KCt;
                    CtKC=CtKCt;
                else
                    mualpha=mualpha/2;
                end
            end          
        end
    end
    if delta~=0
        C=C*diag(alphaC);
    end

  

