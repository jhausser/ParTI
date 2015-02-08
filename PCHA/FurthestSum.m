function i=FurthestSum(K,noc,i,exclude)

% FurthestSum algorithm to efficiently generate initial seeds/archetypes
%
% Usage
%   i=FurthestSum(K,noc,i,exclude)
%
% Input:
%   K           Either a data matrix or a kernel matrix
%   noc         number of candidate archetypes to extract
%   i           inital observation used for to generate the FurthestSum
%   exclude     Entries in K that can not be used as candidates
%
% Output:
%   i           The extracted candidate archetypes

if nargin<4
    exclude=[];
end

[I,J]=size(K);
index=1:J;
index(exclude)=0;
index(i)=0;
ind_t=i;
sum_dist=zeros(1,J);
if J>noc*I % Fast implementation for large scale number of observations. Can be improved by reusing calculations
    Kt=K;    
    Kt2=sum(Kt.^2);
    for k=1:noc+10            
        if k>noc-1 % Remove initial seed        
           Kq=Kt(:,i(1))'*Kt;        
           sum_dist=sum_dist-sqrt(Kt2-2*Kq+Kt2(i(1)));                                    
           index(i(1))=i(1);
           i(1)=[];                       
        end   
        t=find(index);
        Kq=Kt(:,ind_t)'*Kt;        
        sum_dist=sum_dist+sqrt(Kt2-2*Kq+Kt2(ind_t));                            
        [val,ind]=max(sum_dist(t));        
        ind_t=t(ind(1));
        i=[i t(ind(1))];   
        index(t(ind(1)))=0;
    end    
else
    if I~=J || sum(sum(K-K'))~=0 % Generate kernel if K not a kernel matrix
       Kt=K;
       K=Kt'*Kt;        
       K=sqrt(repmat(diag(K)',J,1)-2*K+repmat(diag(K),1,J));            
    end
    Kt2=diag(K)';
    for k=1:noc+10                            
        if k>noc-1
            sum_dist=sum_dist-sqrt(Kt2-2*K(i(1),:)+Kt2(i(1)));
            index(i(1))=i(1);
            i(1)=[];                    
        end           
        t=find(index);
        sum_dist=sum_dist+sqrt(Kt2-2*K(ind_t,:)+Kt2(ind_t));
        [val,ind]=max(sum_dist(t));        
        ind_t=t(ind(1));
        i=[i ind_t];        
        index(t(ind(1)))=0;
    end
end