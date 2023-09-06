function dist = sq_dist(P,Q,Bin_width)
%SQ_DIST Summary of this function goes here
%   Detailed explanation goes here
dist=sum(sum((P-Q).^2))*Bin_width^2;%Integral done with Square Rule assumming uniform sampling
end

