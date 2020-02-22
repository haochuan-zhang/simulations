function  MSE_error=AMP_SE(Input)

%% Load parameters
N=Input.N;
M=Input.M;
alpha=M/N;
nuw=Input.nuw;
IterNum=Input.IterNum;

Gaussian=@(x,m,v) 1./sqrt(2*pi*v).*exp(-(x-m).^2./(2*v));

MSE=1;

for kk=1:IterNum
    
    gamma=1/(nuw+MSE/alpha);
    MSE=1-integral(@(z) tanh(gamma+sqrt(gamma)*z).*Gaussian(z,0,1),-Inf,Inf);
    MSE_error(kk,1)=MSE;
end
end