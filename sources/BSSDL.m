function [se,X,WM,piA,Ae,D] = BSSDL(s,A,n_sample,num_atoms,m,n_iter,OMP_param,SNR_db)

%This function simulate the following process: 

%       1. The receiver gets the mixed signals a.k.a "X" matrix;
%       2. Decomposing the mixed signals "X" into sparse vectors W_m (sparsifying mixed signals);
%       3. Mixing matrix estimation (sparse vectors W_m must be used);
%       4. Separating the sparse decompositions;
%       5. Reconstructing sources from sparse vectors.

%Inputs:
%       s - matrix of num_sources x num_sample.
%       n_sample - sample length in which each column is splitted.
%       num_atoms - Size of dictionary (No. of atoms)
%       optional parameters:
%           m - type of dictionary learning 0: MOD, 1: ODL
%           whether this parameter is greater than 1 the function will use 
%           one of the following "static dictionary": 2: Gabor, 3: DCT,
%           4: Wavelet (level 5)
%           n_iter - number of iteration for MOD

%Outputs:
%       se - estimated sources
%       X - mixes signals alias mixtures
%       WM - sparse representation of the mixture
%       A - input mixing matrix
%       W - demixing matrix
%       D - Dictionary


%Author: ...
%Version: ..
%Date: ...


% track length (every track is supposed to have the same length)
[num_sources,track_length] = size(s);

%% Splitting audio (Blocking stage)
for i = 1:num_sources
    sample(:,:,i) = splitaudio(s(i,:)',n_sample);
end

[patch_length,num_patches] = size(sample(:,:,1));

%% Matrix of unrolled sources
S = zeros(num_sources,patch_length*num_patches);
for i = 1:num_sources
    S(i,:) = reshape(sample(:,:,i),1,[]);
end

%% TODO: create function for creating of the training set
%% Creating training set (different from the original one)
% [c1,~] = audioread('/Users/davidenardone/Dropbox/MATLAB/SM/spams-matlab/BSSproject/input/audio/dataset1/source1.wav');
% [c2,~] = audioread('/Users/davidenardone/Dropbox/MATLAB/SM/spams-matlab/BSSproject/input/audio/dataset1/source2.wav');
% [c3,~] = audioread('/Users/davidenardone/Dropbox/MATLAB/SM/spams-matlab/BSSproject/input/audio/dataset1/source3.wav');
% [c4,~] = audioread('/Users/davidenardone/Dropbox/MATLAB/SM/spams-matlab/BSSproject/input/audio/dataset1/source4.wav');
% t(1,:) = c1(:,1);
% t(2,:) = c2(:,1);
% t(3,:) = c3(:,1);
% t(4,:) = c4(:,1);
% 
% for i = 1:num_sources
%     t_sample(:,:,i) = splitaudio(t(i,:)',n_sample);
% end

train = []; %from the same dataset
% t_train = []; %from a different dataset
for i = 1:num_sources
    train = [train sample(:,:,i)];
%     t_train = [t_train t_sample(:,:,i)];
end

%% Constructing Dictionary 
if m <= 2
    % Learning Dictionary from the sample (in thise case all sample have been taken)
    [D, ~, ~] = process_audio(train,num_sources,num_atoms,m,n_iter);
else % static dictionary
    if m == 3
        [D, ~] = constructDictionaryGabor(n_sample);
    end
    % http://it.mathworks.com/help/wavelet/ref/wmpdictionary.html
    if m == 4
        D = wmpdictionary(n_sample,'lstcpt',{'dct'});
    end
    if m == 5 %Haar wavelet packet
        D = wmpdictionary(n_sample,'lstcpt',{'wphaar'}); % level 5
    end
    if m == 6
        D = wmpdictionary(n_sample,'lstcpt',{'sin','cos'});
    end
    if m == 7
        D = wmpdictionary(n_sample,'lstcpt',{'RnIdent'});
    end
end

% ImD=displayPatches(abs(D));
% subplot(1,3,1);
% imagesc(ImD); colormap('gray');

[~, num_atoms] = size(D);
D = full(D);

% W = zeros(num_atoms,num_patches,num_sources);
% W_s = zeros(num_sources,num_atoms*num_patches);
% ssample = zeros(num_sources,patch_length*num_patches);

% OBSOLETE ...
%% Sparse representation on the Dictionary just learned
% for i = 1:num_sources
%     [w_i, ssample_i] = AudioProcessUsingDict(sample(:,:,i),D,m,OMP_param);
%     W(:,:,i) = w_i; %sparse representation of sample(:,:,i)
% %     W_s(1,:) = reshape(w_i,1,[]);
%     ssample(i,:) = reshape(ssample_i,1,[]);
% end



%%  At this part, the sensor nodes trasmit the signals across the network:
%   ***************************SIMULATION************************

%       1. The receiver gets the mixed signals a.k.a "X" matrix;
%       2. Decomposing the mixed signals "X" into sparse vectors W_m (sparsifying mixed signals);
%       3. Mixing matrix estimation (sparse vectors W_m must be used);
%       4. Separating the sparse decompositions;
%       5. Reconstructing sources from sparse vectors.

%% 1. Mixing process: the receiver gets the mixed signals a.k.a "X" matrix
% A = abs(randn(num_sources,num_sources)); %Random Mixing Vector A
A = abs(A);
A = A./repmat(sqrt(sum(A.^2)),[size(A,1) 1]); %Normalizing Columns of A
X  = mexCalcXY(A,S);

% X = A*S + N, adding noise
if nargin > 7
    X = X + std2(X)*10^(-SNR_db/20)*randn(size(X));
end

%% 2. Decomposing the mixed signals into sparse vectors (Sparsifying the mixtures)
mix = zeros(patch_length,num_patches,num_sources);
WM = zeros(num_sources,num_atoms*num_patches);
% X_est = zeros(num_sources,size(X,2));
for i = 1:num_sources
    mix(:,:,i) = splitaudio(X(i,:),n_sample);
    [wm_i, w_i] = AudioProcessUsingDict(mix(:,:,i),D,OMP_param);
    nnz(wm_i)/(size(wm_i,1)*size(wm_i,2))
    WM(i,:) = reshape(wm_i,1,[]); % concatenating the sparse decompositions of the mixed signals
    
%     X_est(i,:) = reshape(w_i,1,[]); % (unnecessary) reconstructing the mixture signal using the dictionary D
end

% sparsity factor: that is how much sparse is the matrix WM (values close to 0 mean a lot of zero elements)
% sf = nnz(WM)/prod(size(WM));

% dimensionality factor: that how much sparse coefficients we use to
% represent the mixture X (values close to 0 mean high compress factor)
% df = size(WM,2)/size(X,2);

%% 3. Mixing matrix estimation (actually we're using X);
% [~,j,~] = find(WM);
% WM2 = WM(:,unique(j)); % removing columns completely equal to zero

%% method "fgmca"
% [~,j,~] = find(WM);
% WM2 = WM(:,unique(j));
[piA,~] = fgmca(WM,num_sources,500,3);
Ae = abs(inv(piA)); % estimated mixing matrix

% Accuracy estimate mixing matrix:
% https://books.google.it/books?id=hQeaAgAAQBAJ&pg=PT98&lpg=PT98&dq=accuracy+mixing+matrix+blind+source+separation&source=bl&ots=f7rKjz0gH5&sig=SyiWYkX84M-QqHR1WMQXRrC8cVs&hl=it&sa=X&ved=0ahUKEwidwb2x9ozNAhUKuBoKHRkGAHYQ6AEIXDAI#v=onepage&q=accuracy%20mixing%20matrix%20blind%20source%20separation&f=false
% norm(eye(num_sources,num_sources)-(Ae*A))
%% method 1 fuzzy C-means
% [~,j,~] = find(WM);
% WM2 = WM(:,unique(j));
% [W,~] = fcm(full(WM2'),num_sources);
% W = full(W);
% W = W./repmat(sqrt(sum(W.^2)),[size(W,1) 1]);
% W = abs(W);

%% method 2 k-means
% [~,j,~] = find(WM);
% WM2 = WM(:,unique(j));
% [~,W] = kmeans(full(WM2'),num_sources);
% W = full(W);
% W = W./repmat(sqrt(sum(W.^2)),[size(W,1) 1]);
% W = abs(W);
% piA = (W'*W)\W';


%% 4. Separating the sparse decompositions
% param.L = num_sources;
param.eps = 1e-6;
w_rec = zeros(num_sources,size(WM,2));
[~,j,~] = find(WM);
j = unique(j);
for i = 1:numel(j)
    w_rec(:,j(i)) = mexOMP(full(WM(:,j(i))),abs(Ae),param);
end

SS = zeros(num_atoms,num_patches,num_sources);
se = zeros(num_sources,track_length);

%% 5. Reconstruction sources from sparse vectors
for i = 1:num_sources
    SS(:,:,i) = splitaudio(w_rec(i,:),num_atoms);
    tmp = reshape(D*SS(:,:,i),1,[]);
    se(i,:) = tmp(1:track_length);
end

end