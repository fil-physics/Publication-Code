function [p1, dp1, relres1, p0]=ptb_Fit_ppm_complex(M,TE)
%   [p1, dp1, relres, p0]=ptb_Fit_ppm_complex(M)
%    
%   output
%   p1 - field map, may need further unwrapping
%   dp1 - a priori error estimate
%   relres - relative residual
%   p0 - initial phase
%
%   input
%   M - a multi-echo and could be a multi-channel dataset
%       echo needs to be the 4th dimension
%       channel needs to be the 5th dimension
%
%   When using the code, please cite 
%   T. Liu et al. MRM 2013;69(2):467-76
%   B. Kressler et al. IEEE TMI 2010;29(2):273-81
%   de Rochefort et al. MRM 2008;60(4):1003-1009
%
%   The coil combination method is similar to
%   MA. Bernstein et al. MRM 1994;32:330-334
%
%   Adapted from a linear fitting created by Ludovic de Rochefort
%   Modified by Tian Liu on 2011.06.01
%   Modified by Alexey Dimov on 2016.05.12
%   Last modified by Julio Acosta-Cabronero

%Modification to handle one echo datasets - assuming zero phase at TE = 0;
%- AD, 05/12/2016
if size(M,4) == 1
    M = cat(4,abs(M),M);
    single_echo = true;
else
    single_echo = false;
end

if size(M,5)>1
% combine multiple coils together, assuming the coil is the fifth dimension
    M = sum(M.*conj( repmat(M(:,:,:,1,:),[1 1 1 size(M,4) 1])),5);  
    M = sqrt(abs(M)).*exp(1i*angle(M));
end

M= conj(M);
s0=size(M);
L_s0=length(s0);
nechos=size(M,L_s0);

M=reshape(M,[prod(s0(1:L_s0-1)),s0(L_s0)]);
s=size(M);

Y=angle(M(:,1:min(3,nechos)));
% estimate the slope
c=((Y(:,2)-Y(:,1)));
[m ind]=min([abs(c-2*pi),abs(c),abs(c+2*pi)],[],2);
c(ind==1)=c(ind==1)-2*pi;
c(ind==3)=c(ind==3)+2*pi;

for n=1:min(2,nechos-1)
    cd=((Y(:,n+1)-Y(:,n)))-c;
    Y(cd<-pi,(n+1):end)=Y(cd<-pi,n+1:end)+2*pi;
    Y(cd>pi,(n+1):end)=Y(cd>pi,n+1:end)-2*pi;
end
clear c cd % jac
    
A = [1 0;1 1;1 2];
ip = A(1:min(3,nechos),:)\Y(:,1:min(3,nechos))';
p0 = ip(1,:)';
p1 = ip(2,:)';
clear ip Y % jac

dp1 = p1;
tol = norm(p1(:))*1e-4;
iter = 0;
max_iter = 30;

% JAC
% If multi echo, estimate magnitude for SNR weighting
single_echo = true; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JAC temp
if ~single_echo
    B1 = log(M);
    A1 = cat(2,ones(size(TE(1:nechos)')),TE(1:nechos)');
    X1 = A1\B1(:,1:nechos)';
    clear A1 B1
    S0 = exp(X1(1,:)');
    T2star = -1./X1(2,:)';
    clear X1
    for x=1:nechos
        Mest(:,x) = S0.*exp(-TE(x)./T2star);
    end
    clear S0 T2star
else
    Mest = M;
end

% weigthed least square
% calculation of WA'*WA
v1=ones(1,nechos);
v2=(0:(nechos-1));
a11=sum(abs(Mest).^2.*(ones(s(1),1)*(v1.^2)), 2);  % jac
a12=sum(abs(Mest).^2.*(ones(s(1),1)*(v1.*v2)),2);  % jac
a22=sum(abs(Mest).^2.*(ones(s(1),1)*(v2.^2)), 2);  % jac
% inversion
d=a11.*a22-a12.^2;
ai11=a22./d;
ai12=-a12./d;
ai22=a11./d;
clear a11 a12 a22 d % jac

while ((norm(dp1)>tol) &&(iter<max_iter))
    iter = iter+1;
    W = abs(Mest).*exp(1i*(p0*v1 + p1*v2) ); % jac

    % projection
    pr1=sum(conj(1i*W).*(ones(s(1),1)*v1).*(M-W),2);
    pr2=sum(conj(1i*W).*(ones(s(1),1)*v2).*(M-W),2);

    dp0=real(ai11.*pr1+ai12.*pr2);
    dp1=real(ai12.*pr1+ai22.*pr2);
    dp1(isnan(dp1))=0;
    dp0(isnan(dp0))=0;
    
    %update
    p1 = p1+dp1;
    p0 = p0+dp0;
end

% error propagation
dp1=sqrt(ai22);
dp1(isnan(dp1)) = 0;
dp1(isinf(dp1)) = 0;
clear ai* dp0 pr* W % jac

% relative residual
res = M - abs(Mest).*exp(1i*(p0*v1 + p1*v2));      % jac
relres1 = sum(abs(res).^2,2)./sum(abs(Mest).^2,2); % jac
relres1(isnan(relres1)) = 0;
% relres2 = sqrt(sum(angle(res).^2,2));              % jac
% relres2(isnan(relres2)) = 0;                       % jac

p1(p1>pi)=mod(p1(p1>pi)+pi,2*pi)-pi;
p1(p1<-pi)=mod(p1(p1<-pi)+pi,2*pi)-pi;

p0=reshape(p0,s0(1:L_s0-1));
p1=reshape(p1,s0(1:L_s0-1));
dp1=reshape(dp1,s0(1:L_s0-1));
relres1 = reshape(relres1,s0(1:L_s0-1));           % jac
% relres2 = reshape(relres2,s0(1:L_s0-1));           % jac
