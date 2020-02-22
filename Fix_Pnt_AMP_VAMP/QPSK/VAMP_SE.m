function MSE_error=VAMP_SE(obj,Input)

IterNum=Input.IterNum;
MSE=zeros(IterNum,1);
nuw=Input.nuw;
H=obj.H;
HtH=H'*H;
[~,s,~]=svd(HtH);
lambda=diag(s);
Gaussian=@(x,m,v) 1./sqrt(2*pi*v).*exp(-(x-m).^2./(2*v));
gamma_plus=1;

for jj=1:IterNum
   %% back passing   
   hatv_plus=mean(1./(lambda/nuw+gamma_plus));
   gamma_sub=1/hatv_plus-gamma_plus;
   MSE=1-integral(@(u) tanh(gamma_sub+sqrt(gamma_sub).*u).*Gaussian(u,0,1),-Inf,Inf);
   gamma_plus=1/MSE-gamma_sub;
   MSE_error(jj,1)=MSE;
end

end



