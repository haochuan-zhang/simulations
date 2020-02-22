function   obj=MIMO_system(Input)

N=Input.N;
M=Input.M;
rho=Input.rho;
nuw=Input.nuw;
sigma_X=Input.sigma_X;

%% Generate x
x0 = sqrt(sigma_X)*randn(N,1);            % a dense Gaussian vector
pos=rand(N,1) < rho;
x = x0.*pos;  % insert zeros

%% Real-valued Channel
H=randn(M,N)/sqrt(M); 
if(Input.is_GaussAddUnif==1 & Input.is_GaussAddDiscrete==0)
    scFact = 20; % uniform dist in [-scFact*0.5, +scFact*0.5]
    HUniform = sqrt(scFact)*(rand(M,N)-0.5)/sqrt(M);
    normFact = 1/sqrt(1+scFact/12); % uniform distribuiton variance: scFact/12
    H = normFact* (H + HUniform);
end
if(Input.is_GaussAddUnif==0 && Input.is_GaussAddDiscrete==1)% {-3, -1, +1, +3}*scFact
    scFact= 20; % unifromly distributed on discrete numbers {-3, -1, +1, +3}*scFact
    % 2*randi(4)-5 --> {-3, -1, +1, +3}
    HDiscrete = sqrt(scFact)*(2*randi(4, M,N)-5)/sqrt(M);
    normFact = 1/sqrt(1+scFact*5); % uniform distribuiton variance: scFact/12
    H = normFact* (H + HDiscrete);
end
if(Input.is_GaussAddDiscrete==1 && Input.is_GaussAddUnif==1)
    error('Wrong is_GaussAddDiscrete or is_GaussAddUnif! ')
end
%% Noise
w=sqrt(nuw)*randn(M,1);

%% Uncoded system
y=H*x+w;


%% load parameters
obj.x=x;
obj.y=y;
obj.H=H;
obj.pos=pos;

end