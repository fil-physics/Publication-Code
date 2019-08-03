function kf = m1_cosfilt(L,doff)
% 1-D cosine filter
%
% INPUT
% L     Filter length
% doff  Dropoff extent, i.e. # of elements that will be filtered out
% 
% OUTPUT
% kf    Cosine filter

kf                      = ones(L+2,1);
kf(end-doff-1:end,1)    = cos(linspace(0,pi/2,doff+2));
kf(1:doff+2)            = cos(linspace(pi/2,0,doff+2));
kf                      = kf(2:end-1).^2;
