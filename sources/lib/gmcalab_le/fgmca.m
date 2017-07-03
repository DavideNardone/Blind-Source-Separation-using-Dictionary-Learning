function [piA,S] = fgmca(X,NbSources,nmax,Kmin)

%  Sped-Up GMCA
%
%  INPUT - X : the sparse data - channels * samples
%        NbSources : number of sources to estimate
%        nmax : number of iterations (typically 100)
%        Kmin : last k-mad threshold (typically 3)
%
%  OUTPUT - piA : the pseudo inverse of the mixing matrix
%         S   : estimated/thresholded sparse sources
%
%  REQUIRED TOOLBOX :  Wavelab 850 - http://www-stat.stanford.edu/~wavelab/
%  FACULTATIVE TOOLBOXES : MCALab - http://www.greyc.ensicaen.fr/~jfadili/demos/WaveRestore/downloads/mcalab.html
%                   CurveLab - http://www.curvelet.org/
%
%  Version : 7th of August 2007
%          CEA-Saclay / DAPNIA/SEDI-Sap
%          J.Bobin
%
%

[NbChannels,NbSamples] = size(X);

% scaling data to zero mean (centered)
for ll=1:NbChannels
    X(ll,:) = X(ll,:) - mean(X(ll,:))*ones(size(X(ll,:)));   
end


KS = 15;
DKS = (KS-Kmin)/nmax;  %-- Should be changed - Typical values for an appropriate thresholding

% AA = randn(NbChannels,NbSources);
load('/Users/davidenardone/Dropbox/MATLAB/SM/spams-matlab/BSSproject/estimated_mixing_matrix.mat')

SEst_r = zeros(NbSources,NbSamples);

for pp=1:nmax
    piA = inv(AA'*AA)*AA';
%     piA=(AA'*AA)\AA';
    
    SEst_r = piA*X;
    
    SigmaSources = zeros(1,NbSources);
    
    %thresholding step S = lambda*(AA*X^+)
    for ff = 1:NbSources
        
        SEst_r(ff,:) = SEst_r(ff,:).*( abs(SEst_r(ff,:)) > KS*mad(SEst_r(ff,:)) );

        SigmaSources(ff) = std(SEst_r(ff,:));
        
    end
    
    indd = find(SigmaSources > 1e-9);
    
    if ~isempty(indd)
        
        %updating AA, that is the estimated mixing matrix
        %it may be optimazed decomposing this product into 2 single
        %products and using spam library
        AA(:,indd) = X*SEst_r(indd,:)'*inv(SEst_r(indd,:)*SEst_r(indd,:)');        
%         AA(:,indd) = (X*SEst_r(indd,:)')/(SEst_r(indd,:)*SEst_r(indd,:)');

        for ff = 1:length(indd)
            
            AA(:,indd(ff)) = AA(:,indd(ff))/norm(AA(:,indd(ff)));
            
        end
    end
    
    KS = KS - DKS;
    
end

S = SEst_r;