clear all;clear sound;close all;clc;


dirPath='/Users/davidenardone/Dropbox/MATLAB/SM/spams-matlab/BSSproject/input/audio/';
d=dir(dirPath);
field0='y';
field1='Fs';
field2='length';
f = struct(field0,[],field1,[]);
F=[];

option=1;
option_c=1;

switch option
    
    %% RECONSTRUCT AUDIO SIGNAL
    case 1
        
        %% concatenate audio
        
        fprintf(1,'Reading files and concatenating audio files...\n');
        for k=1:3
            if (~strcmp(d(k).name(1),'.'))
                [y,Fs] = audioread(strcat(dirPath,d(k).name));
                F=[F;y];
            end
        end
        
        len_frame=256;
        Y = splitaudio(F,len_frame);
        
        if (option_c==1)
        %% Sparse representation of an audio signal using ODL
        
            %% testing speech
            [y_test,Fs_test]=audioread('/Users/davidenardone/Dropbox/MATLAB/SM/spams-matlab/BSSproject/input/audio/source2.wav');
            Y_TEST = splitaudio(y_test,len_frame);

            %% ODL
            fprintf(1,'Online Dictionary Learning...\n');
            paramdict.iter=100;
            paramdict.batchsize=5;
            paramdict.mode=1;
            paramdict.K=256;
            paramdict.verbose=false;
            paramdict.lambda=10e-6;
            D=mexTrainDL(Y,paramdict);

            %% OMP
            fprintf(1,'Orthogonal Matching Pursuit...\n');
            paramomp.eps=paramdict.lambda;
            paramomp.L=paramdict.K;
            alpha=mexOMP(Y_TEST,D,paramomp);

            fprintf(1,'SPAM...\n');
            Xhat=D*alpha;

            figure;
            plot(Y_TEST(:),'b');
            axis tight
            title('Original signal');
            figure;
            plot(Xhat(:),'b');
            axis tight
            title('Signal reconstruct');

            % MSE
            figure;
            plot((Y_TEST-Xhat).^2,'r')
            axis tight
            title('MSE');

            fprintf(1,'Number of zero in the spars vector alpha: %d\n', nnz(alpha));

            fprintf(1,'MSE: %d\n', sum(mean(abs(Y_TEST-Xhat).^2)));

            fprintf(1,'\npress any key to listen the reconstructed signal...\n');
            pause();
            sound(Xhat(:),Fs_test) 
            
            
        else if option_c==2
        %% Denoising audion signal        
        
            %% adding noise
            sigma=20/255;

            I=Y+(sigma)*randn(size(Y));
            n=10;
            Y2=mexExtractPatches(I,n);
            me=mean(Y2);
            Y2=bsxfun(@minus,Y2,me);


            %% ODL
            fprintf(1,'Online Dictionary Learning...\n');
            paramdict.iter=100;
            paramdict.batchsize=10;
            paramdict.mode=1;
            paramdict.K=256;
            paramdict.verbose=false;
            paramdict.lambda=n*n*(1.1)*sigma*sigma;
            D=mexTrainDL(Y2,paramdict);

            %% OMP
            fprintf(1,'Orthogonal Matching Pursuit...\n');
            paramomp.eps=paramdict.lambda;
            paramomp.L=paramdict.K;
            alpha=mexOMP(Y2,D,paramomp);

            fprintf(1,'SPAM...\n');
            Xhat=D*alpha+ones(size(D,1),1)*me;
            Ihat=mexCombinePatches(Xhat,I,n);
%             Ihat=mergePatch(Xhat,n,n,size(Y,1),size(Y,2));
            
            figure;
            subplot(1,3,1);
            plot(Y(:),'b');
            subplot(1,3,2);
            plot(I(:),'r');
            subplot(1,3,3);
            plot(Ihat(:),'b');
            fprintf(1,'MSE: %d\n', sum(mean(abs(Y-Ihat).^2)));
            
%             sound(Ihat(:),Fs);
            end
            
        end
        
    case 2
%         col=5;
%         row=ceil(length(d)/col);
%         i=1;
        
        fprintf(1,'Reading files...\n');
        for k=1:13
            if (~strcmp(d(k).name(1),'.'))
                [y,Fs] = audioread(strcat(dirPath,d(k).name));
                f.y=y;
                f.Fs=Fs;
                f.length=length(y);
                F=[F;f];
            end
        end
        
        
        M=max([F(:).length]);
        N=size(F,1);
        
        Y = zeros(M,N);
        data = zeros(M,N);
        
        % padding Y matrix
        for i=1:N
            %     Y(1:F(i).length,i)=dct(F(i).y(:,1)); %dct
            Y(1:F(i).length,i)=F(i).y(:,1); %y(:,1) first channel
%             data(1:F(i).length,i)=F(i).y(:,1);
        end
        
        
        % adding noise
        sigma=10e-3;
        I=Y+(sigma)*randn(size(Y));
        
        %% ODL
        fprintf(1,'Online Dictionary Learning...\n');
        paramdict.iter=100;
        paramdict.batchsize=5;
        paramdict.mode=1;
        paramdict.K=256;
        paramdict.verbose=false;
        paramdict.lambda=0.15;
        D=mexTrainDL(Y,paramdict);
            
        
    case 3
        ...
    otherwise
        disp('option wrong');
        return;
end


