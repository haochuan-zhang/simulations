function MSE_error=AMP_Detector(Input,obj)

%% Load Parameters
N=Input.N;
M=Input.M;
nuw=Input.nuw;
IterNum=Input.IterNum;

H=obj.H;
x=obj.x;
y=obj.y;



%% Parametes Initialization
V_old=ones(M,1);
Z_old=y;   

hatx=zeros(N,1);            %m={m_i}, m_i=theta(check_m,check_v)
varx=ones(N,1);             %v={v_i}; v_i=eta(check_m,check_v)

sqrH=abs(H).^2;
sqrHt=sqrH';
Ht=H';

MSE_error=zeros(IterNum,1);
mes=1;

%% Iteration
for ii=1:IterNum
    if ii>1
        mes=Input.mes;
    end
    V=sqrH*varx;
    V_temp1=1./(nuw+V_old);
    Z=H*hatx-(y-Z_old).*V_temp1.*V;
        
    [V,V_old]=damping(V, V_old, mes);
    [Z,Z_old]=damping(Z, Z_old, mes);

    V_temp2=1./(nuw+V);
    check_v=1./(sqrHt*V_temp2);
    check_m=hatx+check_v.*(Ht*((y-Z).*V_temp2));
    
    [hatx,varx]=estimator_x(Input,check_m,check_v);
    MSE=norm(hatx-x)^2/norm(x)^2;  
    MSE_error(ii,1)=MSE;
end
end

function [x,x_old]=damping(x,x_old,mes)
   x=mes*x+(1-mes)*x_old;
   x_old=x;
end

function [hatx, varx]=estimator_x(Input,m,v)

rho=Input.rho;
sigma_X=Input.sigma_X;

Gau=@(x,a,v) 1./sqrt(2*pi*v).*exp(-1./(2*v).*abs(x-m).^2);

C=(rho*Gau(0,m,v+sigma_X))./...
    (rho*Gau(0,m,v+sigma_X)+(1-rho)*Gau(0,m,v));

hatx=C.*(m*sigma_X)./(v+sigma_X);
Ex2=C.*((v*sigma_X)./(v+sigma_X)+abs((m*sigma_X)./(v+sigma_X)).^2);
varx=Ex2-abs(hatx).^2;

end


