function [BW,B,W]=unica(x,e,w,whitening,MaxIter,stop);
% UNICA Algorithm 
%   (Unbiased qNewton algorithm for Independent Component Analysis).
%   (c) Sergio Cruces, Luis Castedo, Andrzej Cichocki.
% 	http://viento.us.es/~sergio     E-mail: sergio@us.es
%   Version: 2.5                    Last update:   27/12/2002.
% 					                (Version: 1.0, 28/07/2001)
%
%======================================================
%
% PURPOSE:    To extract 'e' independent components performing 
%			  unbiased ICA in presence of correlated
%			  Gaussian noise in the mixture. The algorithm 
%			  alternates a qNewton iteration for the estimation 
%			  of the mixing system with a minimum variance 
%			  distorsionless response criteria which 
%			  eliminates from the outputs the interfering 
%			  components and all the noise that is outside 
%			  the extracted signal subspace.
%
% MANDATORY INPUT ARGUMENTS
%   x  Matrix with the data observations  
%		 with dimension (N.sources x Samples).
%      Each row corresponds to a different sensor.
%
% OPTIONAL INPUT ARGUMENTS
%    e        Number of independent componentes to 
%				  extract simultaneously (default e=1).
%	 w        Vector with the weighting coefficients.
%   whitening When is set to 1 enables prewhitening, 
%             when is set to 0 skips this preprocessing.
%   MaxIter   Maximum number of iterations.
%   stop      Sensivity parameter to stop the convergence. 
%
%
% OUTPUT ARGUMENTS
%   BW  		   Composite ICA transformation.
%   B				Semi-orthogonal transformation.
%   W  		   Prewhitening transformation.
%
%=========================================================
%
% Related Bibliography:
%
% [1] S. Cruces, A. Cichocki, L. Castedo, 
%		"Blind Source Extraction in Gaussian Noise", 
%		In proceedings of the 2nd International Workshop on 
%		Independent Component Analysis and Blind Signal 
%		Separation (ICA'2000), pp. 63--68, Helsinki, Finland,
%       June 2000. (ISBN: 951-22-5017-9)
%
% [2] S. Cruces, "An unified view of BSS algorithms",
%		PhD. Thesis (in Spanish), University of Vigo, 
%		Spain, 1999.
%
%=========================================================

eta0=.9;
[n,T]=size(x);
I=eye(e);
if ~exist('e','var'), e=1;end;
if ~exist('w','var'), w=[0 0 1 0 0 0];end;
if ~exist('alpha','var'), alpha=[1 1 1 1 1 1];end;
if ~exist('whitening','var'), whitening=1; end; 
if ~exist('MaxIter','var'), MaxIter=1000; end; 
if ~exist('stop','var'), stop=1e-3; end;
reinic=2;

% PREWHITENING

x=x-mean(x')'*ones(1,T);
Rxx=x*x'/T;
if whitening, 
  fprintf('Prewhitening... ') 
  W=pinv(sqrtm(Rxx));
  fprintf('Done.\n')
else
  Rxx_=pinv(Rxx);	   
  W=eye(n,n);
  disp('Prewhitening is skipped.') 
end 
x=W*x;   % Prewhiten the data.

% Initialization. 

if e==1 B=rand(1,n);B=B/norm(B);  
else B=eye(e,n);
end
A=B';
C11yx=zeros(e,n);C11yy=zeros(e);D2=zeros(e);
C21yx=zeros(e,n);C12yy=zeros(e);D3=zeros(e);
C31yx=zeros(e,n);C13yy=zeros(e);D4=zeros(e);
C41yx=zeros(e,n);C14yy=zeros(e);D5=zeros(e);
C51yx=zeros(e,n);C15yy=zeros(e);D6=zeros(e);

disp('Extraction...')
disp('[Iteration, Index , Convergence]') 
convergence=0;phi=0;

it=0;
while ~convergence 
  it=it+1;
  
  % Outputs
  y=B*x; 
  
  % Statistics.
  
  y_=conj(y);
  if w(2)
      C11yx=(x*y_.')'/T; 
	  C11yy=B*(C11yx)';
 	  v2=real(diag(diag(C11yy)));
	  D2=sign(v2)*diag(diag(abs(v2)).^(alpha(2)-1));
  end
  if w(3)
     C21yx=(x*(y_.*y).')'/T; 
	  C12yy=B*(C21yx)';
     v3=real(diag(diag(C12yy)));
     D3=sign(v3)*diag(diag(abs(v3)).^(alpha(3)-1));     
  end
  if w(4)    
	  C31yx=(x*(y_.*y.*y_).'-2*(x*y_.')*diag(mean((y.*y_).'))-...
        (x*y.')*diag(mean((y_.*y_).')))'/T; 
     C13yy=B*(C31yx)';
     v4=real(diag(diag(C13yy)));   
     D4=sign(v4)*diag(diag(abs(v4)).^(alpha(4)-1));  
  end;
  if w(5)
     C41yx=(x*(y'.^4)-4*x*y'*diag(mean(y.'.^3))-...
        6*x*(y'.^2)*diag(mean(y.'.^2)))'/T;
     C14yy=B*(C41yx)';
     v5=real(diag(diag(C14yy)));
     D5=sign(v5)*diag(diag(abs(v5)).^(alpha(5)-1));     
  end
  if w(6)
	  C51yx=(x*(y'.^5)-5*x*y'*diag(mean(y.'.^5))-...
         10*x*(y'.^2)*diag(mean(y.'.^3))-...
         10*x*(y'.^3)*diag(mean(y.'.^2))+...
         30*x*y'*diag((mean(y.'.^2)).^2))'/T;     
     C15yy=B*(C51yx)';
	  v6=real(diag(diag(C15yy)));
	  D6=sign(v6)*diag(diag(abs(v6)).^(alpha(6)-1));
  end
   
  CxyD=(...
     w(2)*(D2*C11yx)+...
     w(3)*(D3*C21yx)+...
     w(4)*(D4*C31yx)+...
     w(5)*(D5*C41yx)+...
     w(6)*(D6*C51yx))';
 
   
 
  BDelta= B*CxyD-I;
  
  % Learning step size.
  eta=eta0/(1+fix(it/50)/2); 
  mu=min(2*eta0/(1+eta0*(w(2)+w(3)*2+w(4)*3+w(5)*4+w(6)*5)),...
         eta/(1+eta*norm(BDelta+I,1)));
  % QNewton algorithm.

  A=(1-mu)*A+mu*CxyD;            % Update for the mixing system.
  if whitening                   % Uptate for the extraction system.
     B=pinv(A);                  % Data has been Prewhitened.
  else
     B=pinv(A'*Rxx_*A)*A'*Rxx_;     % No prewhitening has been applied.      
  end;
  
  % Convergence?   
  
  phi(1,it)=(sum(sum(abs(BDelta))))/(e*max((e-1),1));   

  delta=phi(1,it);
  convergence =(delta<stop) | (it==MaxIter);
  
  backB=B;
  %if (trace(B*B')>1e3)& (reinic>=0), %original
    if (trace(B*B')>1e9)& (reinic>=0),%modified by Cichocki
      B=rand(1,n);B=B/norm(B);  
     if reinic==0 | e>1, 
        B=eye(size(backB)); 
        convergence=1; 
        break
     else
        reinic=reinic-1;
     end         
  end;
  
  % Draw evolution. 
  if 0
  	for i=1:e;
     subplot(e,1,i);
     plot(y(i,:),'.');
  	end
  	drawnow;
  end
  fprintf('%d  Index: %6.3f  Convergence: %6.5f\r',...
     it,trace(A'*A), phi(1,it),delta);     
  
  %-------- Interrupt window ---------------
% KS
   pause( 1/100 );
   go_next_step = findobj( 'Tag', 'alg_is_run' );
   if isempty(go_next_step)
      fprintf( '\nUser break.\n\n' );
      break;
   end
%-------- Interrupt window ---------------

end;

C_index=phi(1,it);
%fprintf('End.');
%fprintf(' \n');

% Sort outputs according to their variance
[v,ind]= sort(std(y.'));
B=B(ind,:);
BW=B*W;