clear all;close all;clc
root = '/Users/davidenardone/Desktop/BSSproject/';
sources = strcat(root,'sources/');
resources = strcat(root,'resources/');

%% Opening directory
% SiSEC_2016/UND_2016/dev1/
dirDataset = 'audio/';
dataset = 'SiSEC_2015/UND_2015/dev1/female4/';
dirPathIn = 'input/';
dirPathRes = 'results/';
dirPathOutput = 'output/extract/';
dirMixOutput = 'mixes/';

%% Creating path
dirPathIn = strcat(resources,dirPathIn,dirDataset,dataset);
dirPathRes = strcat(resources,dirPathRes,dirDataset,dataset);
dirPathOutput = strcat(resources,dirPathOutput,dirDataset,dataset);
dirMixOutput = strcat(resources,dirMixOutput,dirDataset,dataset);

cd(sources);

%checking whether the input directory exists
if exist(dirPathIn,'dir') == 0
    fprintf(1,'the input directory "%s" doesnt exists!\n',dirPathIn);
    return;
end

d = dir(dirPathIn);

%% Reading files
% supposing that each track has the same length
num_tracks = 4; %number of audio tracks

ind=1;
for k = 1:length(d)
    if( strcmp(d(k).name(1),'.') )
        continue;
    else
        if (ind-1) == num_tracks
            break;
        end
        [channels,fs] = audioread(strcat(dirPathIn,d(k).name));
        s(ind,:) = channels(:,1);
        % http://stats.stackexchange.com/questions/70553/how-to-verify-a-distribution-is-normalized
        ind = ind+1;
    end
end

info = audioinfo(strcat(dirPathIn,d(k).name));

if info.BitsPerSample == 8
    MAXVAL = 2^8-1;
end
if info.BitsPerSample == 16
    MAXVAL = 2^16-1;
end
if info.BitsPerSample == 32
    MAXVAL = 2^32-1;
end

% normalization step [0-1]
% s = mat2gray(s);

[num_sources,track_length] = size(s);

if num_sources ~=  num_tracks
    fprintf(1,'The No. of tracks picked does not correspond with the No. of tracks read\n');
    return;
end

view = 0;
outer_cycle = 1; %10 %changing parameters such as n_sample (DL)
inner_cycle = 1; % to average the resutl %5

MSE_mean = zeros(outer_cycle,num_sources);
SNR_mean = zeros(outer_cycle,num_sources);
PSNR_mean = zeros(outer_cycle,num_sources);

MER_mean = zeros(outer_cycle,num_sources);
SDR_mean = zeros(outer_cycle,num_sources);
SIR_mean = zeros(outer_cycle,num_sources);
SAR_mean = zeros(outer_cycle,num_sources);

CORR_COEF_mean = zeros(outer_cycle,num_sources);
MMC_mean = zeros(outer_cycle,1);

SF_mean = zeros(outer_cycle,1);
DF_mean = zeros(outer_cycle,1);

EL_TIME_mean = zeros(outer_cycle, 1); %elapsed_time_matrix

N_SAMPLE = zeros(outer_cycle,1);
% stop = [1/10 1/100 1/1000 1/10000 1/100000 1/1000000 1/10000000 1/100000000 1/1000000000 1/10000000000]';

%Random Mixing Matrix A
% A = abs(randn(num_sources,num_sources));
load('/Users/davidenardone/Dropbox/MATLAB/SM/spams-matlab/BSSproject/mixing_matrix.mat');

warning('off','MATLAB:MKDIR:DirectoryExists')
for k = 1:outer_cycle
    
    fprintf('Iteration No. %d\n',k);
    MSE = zeros(inner_cycle,num_sources);
    SNR = zeros(inner_cycle,num_sources);
    PSNR = zeros(inner_cycle,num_sources);
    CORR_COEF = zeros(inner_cycle,num_sources);
    
    MER = zeros(inner_cycle,num_sources);
    SDR = zeros(inner_cycle,num_sources);
    SIR = zeros(inner_cycle,num_sources);
    SAR = zeros(inner_cycle,num_sources);
    
    MMC = zeros(inner_cycle,1);
    
    SF = zeros(inner_cycle,1);
    DF = zeros(inner_cycle,1);
    
    elapsed_time_array = zeros(inner_cycle,1);
    
    m = 0;                          %type of Dictionary Learning
    num_atoms = 470;                %Size of dictionary (No. of atoms)
%     n_sample = k*512;
    n_sample = 1024;
    SNR_db = 30;                    % SNR in dB
    N_SAMPLE(k,1) = n_sample;
    n_iter = 10;                   % No. of iterations of MOD algorithm
    
%     OMP_param.L = floor(25.6*k);    %Maximum 0-norm of each column in sparse representation
%     OMP_param.L = n_sample; 
    OMP_param.eps = 1e-6;           %optional, threshold on the squared l2-norm of the residual

    %% PARAMETER OF No. MAX. ITERATION ERICA AND FASTICA
%     MaxIter = 100*k;

    % for averaging the final results
    for t = 1:inner_cycle

        mse = zeros(1,num_sources);
        snr_val = zeros(1,num_sources);
        psnr = zeros(1,num_sources);
        corrCoef = zeros(1,num_sources);

        fprintf('\tinner_cycle %d/%d...',t,inner_cycle);
        alg = 'DL';

        switch alg
            case 'DL'
                %% BSS using Dictionary Learning
                fprintf(1,'Dictionary Learning processing...');
%                 cond = true;
%                 while(cond == true)
                    tic
%                     se,X,WM,piA,Ae,D
                    [se,mix,sparse_mix,~,Ae,D] = BSSDL(s,A,n_sample,num_atoms,m,n_iter,OMP_param);
                    elapsed_time_array(t,:) = toc;
                    [se,idx] = orderSignal(s,se);
                    
                    %% condition block to avoid results involving matrix with singular working precision
%                     P = eye(num_sources);
%                     P = P(:,idx);
%                     piA = inv(Ae*P); % recomputed pseudo inverse
%                     if (rcond(piA) < 1e-15) || (sum(sum(isinf(piA))) > 0)
%                         cond = true;
%                     else
%                         cond = false;
%                     end
%                 end
            case 'FASTICA'
                %% FASTICA
                fprintf('FASTICA processing...');
                tic
                [se, A, Ae] = fastica (s,'only', 'white','verbose','off','epsilon',1e-6);
                elapsed_time_array(t,:) = toc;
                [se,idx] = orderSignal(s,se);
            case 'ERICA'
                fprintf('ERICA processing...');
                tic
                [BWer,Ber,Wer,Ae,se] = erica(s,num_sources,0,1000,1e-6);
                elapsed_time_array(t,:) = toc;
                [se,idx] = orderSignal(s,se);
            otherwise
                disp('Error on choosing the algorithm')
        end
        
        
        %% Computing corrCoef, MSE, PSNR
        for i = 1:num_sources
            corrCoef(1,i) = sum(diag(flipud(corrcoef(s(i,:),se(i,:)))))/2;
            mse(1,i) = sum(mean(abs(s(i,:)-se(i,:)).^2));
            snr_val(1,i) = snr(s(i,:),se(i,:));
            psnr(1,i) = 20*log10(MAXVAL/sqrt(mse(1,i)));
        end
        
        CORR_COEF(t,:) = corrCoef;
        MSE(t,:) = mse;
        SNR(t,:) = snr_val;
        PSNR(t,:) = psnr;

        %% BSS evaluation (higher is the (dB)value, better is the result)
        %IMP: https://www.math.ucdavis.edu/~aberrian/research/voice_separation/
        % https://hal.inria.fr/inria-00544230/PDF/vincent_TASLP06bis.pdf
        %     fprintf('evaluating BSS...\n');
        [sdr,sir,sar,~] = bss_eval_sources(se,s);
        mer = bss_eval_mix(Ae,A);
        
        %% Mixing Matrix Criterion: MMC = || I_n - piAe*P*A ||_1
        % creating permutation matrix
        P = eye(num_sources);
        P = P(:,idx);
        piA = inv(Ae*P); % recomputed pseudo inverse
        MMC(t,:) = norm( eye(num_sources) - piA*A, 1);
        
        if (isinf(MMC(t,:)) == 1)
            MMC(t,:) = 0;           
        end        
        %row-wise
        MER(t,:) = mer';
        SDR(t,:) = sdr';
        SIR(t,:) = sir';
        SAR(t,:) = sar';
        fprintf(1,'done!\n');
        %         fprintf('press any key to continue...');
        %         pause();
        
    end
    CORR_COEF_mean(k,:) = mean(CORR_COEF);
    MSE_mean(k,:) = mean(MSE);
    SNR_mean(k,:) = mean(SNR);
    PSNR_mean(k,:) = mean(PSNR);
    
      
    MER_mean(k,:) = mean(MER);
    SDR_mean(k,:) = mean(SDR);
    SIR_mean(k,:) = mean(SIR);
    MMC_mean(k,:) = mean(MMC);
    SAR_mean(k,:) = mean(SAR);
    
    SF_mean(:,k) = mean(SF);
    DF_mean(:,k) = mean(DF);
    
    EL_TIME_mean(k) = mean(elapsed_time_array);
     
    clear CORR_COEF MSE SNR PSNR MER SDR SIR SAR MMC SF DF elapsed_time_array;
    clear sample;
end

%% plotting graphs
plotGraphs(s,se,sparse_mix);

%% Saving Outputs
saving_audio_output(mix,se,fs,dirMixOutput,dirPathOutput);

%creating results directory whether it doesnt exits
% if exist(dirPathRes,'dir') == 0
%     mkdir(dirPathRes);
% end
% 
% %move into the result directory
% cd(dirPathRes);
% 
% if(strcmp(alg,'DL'))
%     filename = strcat('MOD_new_mbp_workspace','_Dict_learn_mod=',num2str(m),'_num_atoms_Dict=',num2str(num_atoms),'_num_iter=',num2str(n_iter),'.mat');
%     save(filename);
% else
%     %TODO: take mean of all outercycle results
%     filename = strcat('workspace','_algorithm = ',alg,'.mat');
%     save(filename);
% end