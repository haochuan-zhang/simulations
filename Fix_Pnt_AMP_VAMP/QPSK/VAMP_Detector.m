function MSE=VAMP_Detector(obj,Input)

% load parameters
H=obj.H;
y=obj.y;
N=Input.N;
xo=obj.xo;
x=obj.x;
nuw=Input.nuw;
mes=Input.mes;
IterNum=Input.IterNum;

r2=zeros(N,1);
r2_old=r2;
v2_inv=2*ones(N,1);
v2_inv_old=v2_inv;

MSE=zeros(IterNum,1);

for ii=1:IterNum

    
    % Back Passing
    [hatx2,varx2]=estimator2(y,H,nuw,r2,v2_inv, Input.is_diagU ); 
    v1=max(varx2./(1-varx2.*v2_inv),eps);
    r1=v1.*(hatx2./varx2-r2);
    
    %Forward Passing
    [hatx1,varx1]=estimator1(xo,r1,v1, Input.is_diagU ); 
    MSE(ii,1)=norm(hatx1-x).^2/N;
    varx1=max(varx1,5e-13);
    v2_inv=(v1-varx1)./varx1./v1;
    r2=(hatx1.*v1-r1.*varx1)./varx1./v1;
    
    negldx=v2_inv<0;
    v2_inv(negldx)=v2_inv_old(negldx);
    r2(negldx)=r2_old(negldx);
    
    [v2_inv, v2_inv_old]=damping(v2_inv, v2_inv_old, mes);
    [r2, r2_old]=damping(r2, r2_old, mes);
end
end

function [m,v]=estimator1(xo,check_m,check_v, is_diagU)
log_posterior=bsxfun(@times,-1./check_v,abs(bsxfun(@minus,xo,check_m).^2));
log_posterior=bsxfun(@minus,log_posterior,max(log_posterior));  %防止溢出

posterior=exp(log_posterior); 
posterior=bsxfun(@rdivide,posterior,sum(posterior,2));   %得到标准PDF
m=sum(bsxfun(@times,posterior,xo),2);                    %计算PDF的均值
v=(sum(posterior.*abs(bsxfun(@minus,m,xo).^2),2));       %计算PDF的方差
if (is_diagU==1)
    v=mean(v);
end
end

function [hatx2,varx2]=estimator2(y,A,nuw,r2,v2_inv, is_diagU)
Q=((A'*A)+nuw*diag(v2_inv))^(-1);
hatx2=Q*(nuw*r2+A'*y);
varx2=nuw*real(diag(Q));
if (is_diagU==1)
    varx2=mean(varx2);
end
end

function [x,x_old]=damping(x, x_old, mes)
x=mes*x+(1-mes)*x_old;
x_old=x;
end