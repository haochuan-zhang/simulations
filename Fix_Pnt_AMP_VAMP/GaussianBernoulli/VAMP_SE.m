function MSE_error=VAMP_SE(obj,Input)

IterNum=Input.IterNum;
sigma_X=Input.sigma_X;
rho=Input.rho;
Gaussian=@(x,m,v) 1./sqrt(2*pi*v).*exp(-(x-m).^2./(2*v));
Pm=@(m,v) rho*Gaussian(0,m,v+sigma_X)+(1-rho)*Gaussian(0,m,v)+eps;
nuw=Input.nuw;

MSE=zeros(IterNum,1);
MSE_error=MSE;

H=obj.H;
HtH=H'*H;
[~,s,~]=svd(HtH);
lambda=diag(s);
gamma_plus=1;

for jj=1:IterNum
%    hatv_plus=mean(1./(lambda/nuw+gamma_plus));
%    gamma_sub=1/hatv_plus-gamma_plus;
%    MSE=1-integral(@(u) tanh(gamma_sub+sqrt(gamma_sub).*u).*Gaussian(u,0,1),-Inf,Inf);
%    gamma_plus=1/MSE-gamma_sub;
%    MSE_error(jj,1)=MSE;

hatv_plus=mean(1./(lambda/nuw+gamma_plus));
gamma_sub=1/hatv_plus-gamma_plus;

v=1/gamma_sub;
Ex2Mean=rho*sigma_X;
sqr_ExMean=(rho*sigma_X/(v+sigma_X))^2*...
    integral(@(m) m.^2.*(Gaussian(0,m,v+sigma_X)).^2./Pm(m,v),-Inf,Inf);
MSE=Ex2Mean-sqr_ExMean;

gamma_plus=1/MSE-gamma_sub;
MSE_error(jj,1)=MSE;

end



