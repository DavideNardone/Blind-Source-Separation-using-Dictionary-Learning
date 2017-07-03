function [D,w,y2] = process_audio(audio_data,num_sources,K,mod,noIt)


if nargin<2
    mod=0;
end

x = double(audio_data);
[dim_sigs,n_split] = size(x);

%% Method of Optimal Directions (MOD)
% https://www.amazon.co.uk/dp/144197010X/ref=wl_it_dp_o_pC_S_ttl?_encoding=UTF8&colid=B76SBGBMQBCT&coliid=I2DIS3KTRJBH7V#reader_144197010X

% Task: Train a dictionary A to sparsely represent the data x_i, by
% approximating the solution.
if (mod==0)
    
    % Initialize Dictionary using K randomly columns of x
    r = randi(size(x,2),K,1);
    D = x(:,r);
    D = mexNormalize(D); % Normalize the columns of D

    %Maximum 0-norm of each column in sparse representation
    param.L = 10;     %(not more than L non-zeros coefficients in w)
    param.eps = 1e-6;
    
%     eps = 1e-3;
%     w = mexOMP(x,D,param);
    for it=1:noIt
%         old_rule = norm( (x-D*w),'fro')^2;
        %% Sparse Coding Stage: Use OMP to approximate the solution of:
        % min_{A} ||X-DA||_2^2  s.t. ||A||_0 <= L
        w = mexOMP(x,D,param);
        %% MOD Dictionary Update Stage: Update the dictionary by the formula:
        % min_{D} = ||X-DA||_F^2 = XA'(AA')^-1 = XA^+
        xwt = mexCalcXAt(x,w);
        wwt = mexCalcAAt(w);
        ww_inv = mexInvSym(wwt);
        D = mexCalcXY(xwt,ww_inv);
        D = mexNormalize(D);  %Make each column have unit norm
        new_rule = norm( (x-D*w),'fro')^2;
        fprintf(1,'iteration No. %d: ||x-DA||_F^2 = %f\n', it, new_rule);
        %% Stopping Rule
%         if( new_rule+eps < old_rule )
%             continue;
%         end
    end
    
elseif (mod==1)
%         r = randi(size(x,2),K,1);
%         K=256;
%         D = dictnormalize(x(:,500+(1:K)));
        r = randi(size(x,2),K,1);
        D = x(:,r);
        D = mexNormalize(D); % Normalize the columns of D
%         D = mexNormalize(x(:,randperm(K)));
%         D = x(:,r);
%         D = dictnormalize(D);
        nofIt = 10;
        dlsMOD = struct('D',D, 'Met','K-SVD', 'vsMet','ORMP', 'vsArg',struct('tnz',4));
        dlsMOD.snr = zeros(1,nofIt);
        tic;
        for i=1:nofIt;
            dlsMOD = dlfun(dlsMOD, x, 1); 
        end; 
        toc;
        D = dlsMOD.D;
elseif (mod==2) % Online Dictionary Learning (ODL)
    
    param.iter = noIt;
%     param.batchsize = 400;    %(512 default)
    param.mode = 5;             %(2 default)
    param.K = K;
    param.verbose = false; 
    param.lambda = 0.15;
    D = mexTrainDL(x,param);
    D = mexNormalize(D);  %Make each column have unit norm
end
w=[];
y2=[];
%un-comment later
% w=mexOMP(x,D,param);
% y=D*w;
% y2=reshape(y,1,[]);


end

