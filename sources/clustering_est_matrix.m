clear all;close all;clc;


s = 3; % number of sound
p = 2; % number of micros
n = 16000*5; % number of sample

n_sample = [1,n];

[S(1,:),fs1] = audioread('/Users/davidenardone/Dropbox/MATLAB/SM/spams-matlab/BSSproject/input/audio/SiSEC_2016/UND_2016/dev1/female3/dev1_female3_src_1.wav',n_sample,'double'); %utterance
[S(2,:),fs2] = audioread('/Users/davidenardone/Dropbox/MATLAB/SM/spams-matlab/BSSproject/input/audio/SiSEC_2016/UND_2016/dev2/dev2_nodrums/dev2_nodrums_inst_src_3.wav',n_sample,'double'); %guitar
[S(3,:),fs3] = audioread('/Users/davidenardone/Dropbox/MATLAB/SM/spams-matlab/BSSproject/input/audio/SiSEC_2016/UND_2016/dev2/dev2_wdrums/dev2_wdrums_inst_src_2.wav',n_sample,'double'); %drum


%% Mixing Matrix

A = [0.6118 0.9648 0.2360;0.7910 0.2629 0.9718];

Xn=A*S;


% normalize btw [-1 1]
% for i=1:2
%     
%     maxS = max(X(i,:));
%     minS = min(X(i,:));
%     vx = X(i,:)*2/(max(X(i,:))-min(X(i,:)));
%     
%     Xn(i,:) = vx-(min(vx)+1);   % Normalised vector
% 
% 
% end


% # to "de-normalize", apply the calculations in reverse
% SD = (Sn./2+0.5) * (maxS-minS) + minS;


%% computing STFT
% Xs(:,:,1) = spectrogram(Xn(1,:));
% Xs(:,:,2) = spectrogram(Xn(2,:));



options.n = n;
w = 128;   % size of the window
q = w/4;    % overlap of the window
%% Compute the STFT of the micros, and store them into a matrix |Y|.

for i=1:p
    Xs(:,:,i) = perform_stft(Xn(i,:),w,q, options);
end


mf = size(Xs,1);
mt = size(Xs,2);
P = reshape(Xs, [mt*mf p]);
P = [real(P);imag(P)];


% number of displayed points
npts = 6000; 
% display original points
sel = randperm(n); sel = sel(1:npts);
figure;
plot(Xn(1,sel), Xn(2,sel), '.');
% axis([-1 1 -1 1]*5);
% set_graphic_sizes([], 20);
title('Time domain');

%EXO
%% Display some points of |P| in the transformed (time/frequency) domain.
sel = randperm(size(P,1)); sel = sel(1:npts);
figure;
plot(P(sel,1), P(sel,2), '.');
% axis([-1 1 -1 1]*5);
% set_graphic_sizes([], 20);
title('Transformed domain');
%EXO



%normalize
X1n = Xs(:,:,1)/norm(Xs(:,:,1));
X2n = Xs(:,:,2)/norm(Xs(:,:,2));

x = [abs(X1n(:)) , abs(X2n(:))];

% p1 = Xn(:,1:n_sample(2)/3);
% p2 = Xn(:,n_sample(2)/3+1:2*n_sample(2)/3);
% p3 = Xn(:,2*n_sample(2)/3+1:3*n_sample(2)/3);
% 
% mp1=mean(mean(p1));
% mp2=mean(mean(p2));
% mp3=mean(mean(p3));
% M=[mp1;mp2;mp3];


% [~,C] = kmeans(x,3,'Start',[M M);
[~,C] = kmeans(x,3);

A
C'