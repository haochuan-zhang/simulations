function   obj=MIMO_system(Input)

mod_size=Input.mod_size;
N=Input.N;
M=Input.M;
nuw=Input.nuw;

%% Generate x
sym = modem.qammod(2^mod_size);
normal_scal = 1/sqrt((2/3)*(2^mod_size-1)) ;      %QAM normalization 
sym.input='bit';                                  %the type of input data shoud be 'bit'
sym.symbolorder='gray';
informationBit = round(rand(N*mod_size,1)) ;      %产生串行二进制序列
informationSym = modulate(sym, informationBit);   %对二进制序列进行QAM调制
x = normal_scal * reshape(informationSym,N,1) ;   %得到归一化的QAM基带信号

%% Channel
H=(randn(M,N)+1j*randn(M,N))/sqrt(2*M);
if(Input.is_GaussAddUnif==1 && Input.is_GaussAddDiscrete==0)
    scFact = 20; % uniform dist in [-scFact*0.5, +scFact*0.5]
    HUniform = sqrt(scFact)*((rand(M,N)-0.5) +1j*(rand(M,N)-0.5))/sqrt(2*M);
    normFact = 1/sqrt(1+scFact/12); % uniform distribuiton variance: scFact/12
    H = normFact* (H + HUniform);
end
if(Input.is_GaussAddUnif==0 && Input.is_GaussAddDiscrete==1)% {-3, -1, +1, +3}*scFact
    scFact= 20; % unifromly distributed on discrete numbers {-3, -1, +1, +3}*scFact
    % 2*randi(4)-5 --> {-3, -1, +1, +3}
    HDiscrete = sqrt(scFact)*((2*randi(4, M,N)-5) +1j*(2*randi(4, M,N)-5))/sqrt(2*M);
    normFact = 1/sqrt(1+scFact*5); % uniform distribuiton variance: scFact/12
    H = normFact* (H + HDiscrete);
end
if(Input.is_GaussAddDiscrete==1 && Input.is_GaussAddUnif==1)
    error('Wrong is_GaussAddDiscrete or is_GaussAddUnif! ')
end
%% Noise
w=sqrt(nuw/2)*(randn(M,1)+1j*randn(M,1));   %产生高斯噪声

%% Uncoded system
y=H*x+w;

%% load parameters
obj.x=x;
obj.y=y;
obj.H=H;
obj.xo=Gen_Constellation(mod_size); % Generate original constellation sets
end