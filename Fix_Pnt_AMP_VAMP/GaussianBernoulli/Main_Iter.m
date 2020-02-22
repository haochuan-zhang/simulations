clc;
clear all;
close all

%% Parameters Setting
N=1024;                 %Dimension of x
M=512;                  %dimension of y
mes=1;               %damping factor
Iter_Num=1e2;           %Iteration numbers
IterNum=10;
rho=0.05;                %sparse factor
snr=10;
is_GaussAddUnif=0 % 1: i.i.d. Gauss+Uniform[0~1] H 
                  % 0: i.i.d. Gaussian H.
is_GaussAddDiscrete=0 % 1: i.i.d. Gauss+Discrete[-3,-1,+1,+3] H
                     % 0: i.i.d. Gaussian H.
is_diagU=1; % 1: use diagU (uniform diagonals), 0: use diagV (non-uniform diagonals) of the covariance matrix

%% Load parameters
Input.N=N;
Input.M=M;
Input.mes=mes;
Input.IterNum=IterNum;
Input.nuw=10^(-snr/10);
Input.rho=rho;
Input.sigma_X=1/rho;

Input.is_diagU=is_diagU;
Input.is_GaussAddUnif=is_GaussAddUnif;
Input.is_GaussAddDiscrete=is_GaussAddDiscrete;
%% Array setting
for kk=1:Iter_Num
    
    obj=MIMO_system(Input);
    VAMP_MSE(:,kk)=VAMP_Detector(obj,Input);
    SE_VAMP(:,kk)=VAMP_SE(obj,Input);   
    AMP_MSE(:,kk)=AMP_Detector(Input,obj);   
    if (mod(kk,Iter_Num/10)==0)
        disp(kk/Iter_Num*10);
    end
end

SE_AMP=AMP_SE(Input);
for ii=1:IterNum
    VAMP_MSE_mean(ii,1)=mean(VAMP_MSE(ii,:));
    AMP_MSE_mean(ii,1)=mean(AMP_MSE(ii,:));
    SE_VAMP_mean(ii,1)=mean(SE_VAMP(ii,:));
    SE_AMP_mean(ii,1)=mean(SE_AMP(ii,:));
end

iter=1:IterNum;
semilogy(iter,  VAMP_MSE_mean, 'LineStyle', '-','LineWidth', 1,  'Color','b', 'Marker', 's', 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b' );   
hold on;
semilogy(iter,  AMP_MSE_mean, 'LineStyle', '-','LineWidth', 1,  'Color','b', 'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b' );   
hold on;
semilogy(iter,  SE_VAMP_mean, 'LineStyle', 'none','LineWidth', 1,  'Color','r', 'Marker', 'x', 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'r' );   
hold on;
semilogy(iter,  SE_AMP_mean, 'LineStyle', 'none','LineWidth', 1,  'Color','r', 'Marker', '+', 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'r' );   
hold on;

legend('VAMP-Algo','AMP-Algo','VAMP-SE', 'AMP-SE'); hold on;
xlabel('iteration');
ylabel('MSE(dB)');
saveas(figure(1), ['H_Unif',num2str(is_GaussAddUnif),'_Disrc',num2str(is_GaussAddDiscrete),'_',num2str(M),'x',num2str(N),'m',num2str(mes),'_',num2str(randi(1e6)),'.fig'])