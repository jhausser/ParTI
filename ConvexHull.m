function [ vol ] = ConvexHull( data )
%A small aid function that calculates the convex hull of data but also
% handles degenerate cases.
dim  = size(data,2);
switch dim 
    case 0
        vol =  1E-10;
    case 1 
        vol = max(data)-min(data);
    otherwise
        [~, vol] = convhulln(data);
end

