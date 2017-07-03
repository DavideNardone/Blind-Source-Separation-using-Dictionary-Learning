% --- Bi-Orthogonal wavelet transform (Wavelab is mandatory)
% --- Denoised Sources are direct outputs of the blind-GMCA algorithm
% --- Raw Sources are obtained by applying the pseudo-inverse of the estimated mixing matrix to the data

clear all;close all

SNR_db = 30;  %-- SNR in dB

nc = 4;  %-- Number of channels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im1 = double(imread('paraty256.tif'));
im2 = double(imread('paraty256_2.tif'));
im3 = double(imread('sosdelrey.tif'));
im4 = double(imread('pakhawaj.tif'));

im1 = im1 - mean(reshape(im1,1,numel(im1)))*ones(size(im1));
im2 = im2 - mean(reshape(im2,1,numel(im2)))*ones(size(im2));
im3 = im3 - mean(reshape(im3,1,numel(im1)))*ones(size(im1));
im4 = im4 - mean(reshape(im4,1,numel(im2)))*ones(size(im2));


%wavelet coefficients
qmf = MakeONFilter('Battle',5);
wc1 = FWT2_PO(double(im1),1,qmf);
wc2 = FWT2_PO(double(im2),1,qmf);
wc3 = FWT2_PO(double(im3),1,qmf);
wc4 = FWT2_PO(double(im4),1,qmf);

% overdetermined case n>m (i.e. 10>4)
A = randn(nc,4);

%--- input data for fgmca
Xw = A*[reshape(wc1,1,numel(wc1));reshape(wc2,1,numel(wc2));reshape(wc3,1,numel(wc3));reshape(wc4,1,numel(wc4))]; 

Xw = Xw + std2(Xw)*10^(-SNR_db/20)*randn(size(Xw));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[piA,S] = fgmca(Xw,4,100,3);

Sw =  piA*Xw;

wwc1 = reshape(Xw(1,:),length(im1),length(im1));
wwc2 = reshape(Xw(2,:),length(im1),length(im1));
wwc3 = reshape(Xw(3,:),length(im1),length(im1));
wwc4 = reshape(Xw(4,:),length(im1),length(im1));

X1 = IWT2_PO(wwc1,1,qmf);
X2 = IWT2_PO(wwc2,1,qmf);
X3 = IWT2_PO(wwc3,1,qmf);
X4 = IWT2_PO(wwc4,1,qmf);

figure
subplot(221)
imnb(X1)
title('Mixture 1')
subplot(222)
imnb(X2)
title('Mixture 2')
subplot(223)
imnb(X3)
title('Mixture 3')
subplot(224)
imnb(X4)
title('Mixture 4')

wwc1 = reshape(S(1,:),length(im1),length(im1));
wwc2 = reshape(S(2,:),length(im1),length(im1));
wwc3 = reshape(S(3,:),length(im1),length(im1));
wwc4 = reshape(S(4,:),length(im1),length(im1));

Sdn1 = IWT2_PO(wwc1,1,qmf);
Sdn2 = IWT2_PO(wwc2,1,qmf);
Sdn3 = IWT2_PO(wwc3,1,qmf);
Sdn4 = IWT2_PO(wwc4,1,qmf);

wwc1 = reshape(Sw(1,:),length(im1),length(im1));
wwc2 = reshape(Sw(2,:),length(im1),length(im1));
wwc3 = reshape(Sw(3,:),length(im1),length(im1));
wwc4 = reshape(Sw(4,:),length(im1),length(im1));

S1 = IWT2_PO(wwc1,1,qmf);
S2 = IWT2_PO(wwc2,1,qmf);
S3 = IWT2_PO(wwc3,1,qmf);
S4 = IWT2_PO(wwc4,1,qmf);


figure
subplot(221)
imnb(Sdn1)
title('Thresholded Source - GMCA Output - 1')
subplot(222)
imnb(Sdn2)
title('Thresholded Source - GMCA Output - 2')
subplot(223)
imnb(Sdn3)
title('Thresholded Source - GMCA Output - 3')
subplot(224)
imnb(Sdn4)
title('Thresholded Source - GMCA Output - 4')

figure
subplot(221)
imnb(S1)
title('Raw Source - 1')
subplot(222)
imnb(S2)
title('Raw Source - 2')
subplot(223)
imnb(S3)
title('Raw Source - 3')
subplot(224)
imnb(S4)
title('Raw Source - 4')
