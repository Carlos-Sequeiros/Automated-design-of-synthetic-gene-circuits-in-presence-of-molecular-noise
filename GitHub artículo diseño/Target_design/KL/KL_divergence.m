function KL = KL_divergence(P,Q,Bin_width)
%KL_DIVERGENCE Summary of this function goes here
%   Detailed explanation goes here
f=P.*log(P./Q);%Function to integrate
f(isnan(f))=0;%solve limit indeterminations
f(isinf(f))=0;
KL=sum(sum(f))*Bin_width^2;%Integral done with Square Rule assumming uniform sampling
end

