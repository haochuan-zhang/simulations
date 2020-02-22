function MSE_error=AMP_Detector(Input,obj)

%% Load Parameters
N=Input.N;
M=Input.M;
nuw=Input.nuw;
IterNum=Input.IterNum;
mes=Input.mes;

H=obj.H;
x=obj.x;
xo=obj.xo;
y=obj.y;



%% Parametes Initialization
V_old=ones(M,1);
Z_old=zeros(M,1);   

m=zeros(N,1);            %m={m_i}, m_i=theta(check_m,check_v)
v=ones(N,1);             %v={v_i}; v_i=eta(check_m,check_v)

sqrH=abs(H).^2;
sqrHt=sqrH';
Ht=H';

MSE_error=zeros(IterNum,1);
MSE_old=1;

%% Iteration
for ii=1:IterNum
    
    V=sqrH*v;
    V_temp1=1./(nuw+V_old);
    Z=H*m-(y-Z_old).*V_temp1.*V;
        
    [V,V_old]=damping(V, V_old, mes);
    [Z,Z_old]=damping(Z, Z_old, mes);

    V_temp2=1./(nuw+V);
    check_v=1./(sqrHt*V_temp2);
    check_m=m+check_v.*(Ht*((y-Z).*V_temp2));
    
    [m,v]=estimator_x(xo,check_m,check_v);
    MSE=norm(m-x).^2/N;
    if MSE>MSE_old
        MSE_error(ii:IterNum,1)=MSE;
        break;
    end
    MSE_old=MSE;
    MSE_error(ii,1)=MSE;
end
end

function [x,x_old]=damping(x,x_old,mes)
   x=mes*x+(1-mes)*x_old;
   x_old=x;
end

function [umean,uvar]=estimator_x(xo,v,wvar)
logpxr = bsxfun(@times, -1./wvar, abs(bsxfun(@minus, v, xo)).^2);
logpxr = bsxfun(@minus, logpxr, max(logpxr) );            
pxr = exp(logpxr);
pxr = bsxfun(@rdivide, pxr, sum(pxr,2) );
umean = sum(bsxfun(@times, pxr, xo), 2); 
uvar = sum(pxr .* abs(bsxfun(@minus, umean, xo)).^2, 2);
end


