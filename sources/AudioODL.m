clear all;close all;clc;
M=1000; %number of data point
N=10; %number of sample
x=linspace(-5,5,M); %function domain
A=1; %amplitude

noiseAmplitude = .1;
j=1;
noisy_y=0;

% y=A*sin(2.*x)+cos(3.*x);

% [y,Fs] = audioread('/Users/davidenardone/Dropbox/MATLAB/SM/dataset/AUDIO/PIANO.WAV');
% y=y';

Y=zeros(length(x),N);



for i=1:N
    
    y=A*cos(i*x);
    noisy_y = (y + noiseAmplitude * randn(1, length(y)));
    Y(:,i) = noisy_y;
    noiseAmplitude = noiseAmplitude + 1e-3;
    
    subplottight(10,10,i), plot(1:length(y), Y(:,i));
%     plot(x,noisy_y);
    j = j+1;
    
end

y_test=sin(7*x);
Y_TEST=repmat(y_test',1,N);

paramdict.iter=100;
paramdict.batchsize=100;
paramdict.mode=1;
paramdict.K=256;
paramdict.verbose=false;
% paramdict.lambda=n*n*(1.1)*sigma*sigma;
paramdict.lambda=10e-6;
D=mexTrainDL(Y,paramdict);

paramomp.eps=paramdict.lambda;
paramomp.L=paramdict.K;
alpha=mexOMP(Y_TEST,D,paramomp);

Xhat=D*alpha;
% Ihat=mexCombinePatches(Xhat,I,n);

figure;
% for i=1:N
%     subplottight(10,10,i), plot(1:size(Xhat,1), Xhat(:,i)); 
% end
plot(x,Xhat);

% nnz(alpha)
% mean(abs(noisy_y-Xhat).^2)
% figure;
% plot((noisy_y-Xhat).^2,'r')
