function d = distFromPoly(P,arch)

% arch is a matrix with archtypes as columns
% P is a point, each point is a row
% P0,arch,Alpha,dim,Points,PointsT
dim= size(arch,2)-1;
archR=arch(:,2:end);


if dim==0
    d=sum((P'-arch(:,1)).^2).^0.5;
else
    archT=bsxfun(@minus,archR,arch(:,1));
    inProd=archT'*archT;
    P0=P'-arch(:,1);
    Alpha=(P0'*archT)*inv(inProd);
    Ptag=archT*Alpha';
    ind=Alpha>0;
    if (sum(Alpha)<=1 & all(Alpha>=0))
        d=sum((Ptag-P0).^2).^0.5;
    elseif any(Alpha<0)
        d=distFromPoly(P,[arch(:,1),archR(:,ind)]);
    else
        d=distFromPoly(P,archR);
    end        
end
    