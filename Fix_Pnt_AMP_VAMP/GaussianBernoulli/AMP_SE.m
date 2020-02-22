function  MSE_error=AMP_SE(Input)

%% Load parameters
N=Input.N;
M=Input.M;
alpha=M/N;
nuw=Input.nuw;
IterNum=Input.IterNum;
sigma_X=Input.sigma_X;

Gaussian=@(x,m,v) 1./sqrt(2*pi*v).*exp(-(x-m).^2./(2*v)); %ÊµGaussian·Ö²¼

MSE=1;
rho=Input.rho;

Pm=@(m,v) rho*Gaussian(0,m,v+sigma_X)+(1-rho)*Gaussian(0,m,v)+eps;

for kk=1:IterNum
    v=nuw+1/alpha*MSE;
    Ex2Mean=rho*sigma_X;
    sqr_ExMean=(rho*sigma_X/(v+sigma_X))^2*...
        integral(@(m) m.^2.*(Gaussian(0,m,v+sigma_X)).^2./Pm(m,v),-Inf,Inf);
    MSE=Ex2Mean-sqr_ExMean;
    MSE_error(kk,1)=MSE;
end
end