function [B,W,C_index]=simbec(x,e,w,alpha,MaxIter,whitening,stop);
%  SIMBEC Algorithm 
%  (Simultaneous Blind signal Extraction using Cumulants).
% 	(c) Sergio Cruces, Andrzej Cichocki, S-i. Amari.
%     E-mail: sergio@us.es
%   Version 2.0                     Last update   27/12/2002.
% 	                                Version: 1.0, 28/07/2001.
%
%======================================================
%
% PURPOSE: Perform ICA using simultaneous extraction 
%			  of 'e' independent components that are present
%			  in the given data.
%
% OUTPUT ARGUMENTS
%   B  		    Unitary ICA transformation.
%   W  		    Prewhitening transformation.
%   C_index     Final value of the contrast function.
%
% MANDATORY INPUT ARGUMENTS
%   x  Matrix with the data observations  
%		 with dimension (N.sources x Samples).
%      Each row corresponds to a different sensor.
%
% OPTIONAL INPUT ARGUMENTS
%   e         Number of independent componentes to 
%				  extract simultaneously (default e=1).
%   w         Vector with the weighting coefficients.
%   alpha     Vector with the exponents for cumulants.
%   MaxIter   Maximum number of iterations.
%   stop      Sensivity parameter to stop the convergence. 
%
%=========================================================
%
% Related Bibliography:
%
% [1] S. Cruces, A. Cichocki, S-i. Amari, "On a new blind 
%       signal extraction algorithm: different criteria 
%       and stability analysis", IEEE Signal Processing 
%       Letters, vol. 9(8), pp. 233 - 236, Aug. 2002.
% [2] S. Cruces, A. Cichocki, S-i. Amari, 
% 		"Criteria for the Simultaneous Blind Extraction of
%		Arbitrary Groups of Sources" in the 3rd international 
%		conference on Independent Component Analysis and 
%		Blind Signal Separation. SanDiego, California, USA, 2001.
% [3] S. Cruces, A. Cichocki, S-i. Amari, 
%	 	"The Minimum Entropy and Cumulant Based Contrast 
%       Functions for Blind Source Extraction",  in Lecture 
%       Notes in Computer Science, Springer-Verlag. 
%       IWANN'2001, vol. II, pp. 786--793, (ISBN: 3-540-42237-4) 
%
% [4] S-i. Amari, "Natural Gradient Learning for over- and 
% 		under-complete bases in ICA", Neural Computation, 
%		Vol. 11, pp. 1875-1883, 1999.
%
%=========================================================

verboseKS = 1;

[n,T]=size(x);
if ~exist('e','var'), e=1;end;
if ~exist('w','var'), w=[0 0 1 0 0 0];end;
if ~exist('alpha','var'), alpha=[1 1 1 1 1 1];end;
if ~exist('whitening','var'), whitening=1; end;
if ~exist('MaxIter','var'), MaxIter=1000; end;
if ~exist('stop','var'), stop=1e-5; end;
reinic=2;
eta0= 1;
I=eye(e);

% PREWHITENING

x=x-mean(x')'*ones(1,T);
if whitening, 
   if verboseKS,
      fprintf('Prewhitening...') 
   end
   
  Rxx=x*x'/T;
  [u,d,v]=svd(Rxx+eye(n));
  d=diag(d)-1; 
  n=max(find(d>1e-14)); 
  W=(u(:,1:n)*diag(real(sqrt(1./d(1:n))))*u(:,1:n)'); 
   if verboseKS,
      fprintf('Done.\n')
   end
else
  W=eye(n,n);
   if verboseKS,
      disp('Prewhitening is skipped.') 
   end
end 
x=W*x;   % Prewhiten the data.

% Initialization. 

if e==1 B=rand(1,n);B=B/norm(B);  
else B=eye(e,n);
end
C11yx=zeros(e,n);C11yy=zeros(e);D2=zeros(e);
C21yx=zeros(e,n);C12yy=zeros(e);D3=zeros(e);
C31yx=zeros(e,n);C13yy=zeros(e);D4=zeros(e);
C41yx=zeros(e,n);C14yy=zeros(e);D5=zeros(e);
C51yx=zeros(e,n);C15yy=zeros(e);D6=zeros(e);

if verboseKS,
	disp('Extraction...')
   disp('[Iteration, Trace(|Cumulant|), convergence]') 
end
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
  
  % Learning step size.
  eta=eta0/(1+fix(it/50)/2); 
  mu=eta/2/max(...
	 w(2)*abs(diag(C11yy).^alpha(2))+...
     w(3)*abs(diag(C12yy).^alpha(3))+...
     w(4)*abs(diag(C13yy).^alpha(4))+...
     w(5)*abs(diag(C14yy).^alpha(5))+...
     w(6)*abs(diag(C15yy).^alpha(6)));   
  
  % Gradient algorithm on the Stifel manifold.
if 0
  B=B+mu*(...
     w(2)*(D2*C11yx-C11yy*D2*B)+...
     w(3)*(D3*C21yx-C12yy*D3*B)+...
     w(4)*(D4*C31yx-C13yy*D4*B)+...
     w(5)*(D5*C41yx-C14yy*D5*B)+...
     w(6)*(D6*C51yx-C15yy*D6*B));
else  
 SS=(...
     w(2)*(D2*C11yx-C11yy*D2*B)+...
     w(3)*(D3*C21yx-C12yy*D3*B)+...
     w(4)*(D4*C31yx-C13yy*D4*B)+...
     w(5)*(D5*C41yx-C14yy*D5*B)+...
     w(6)*(D6*C51yx-C15yy*D6*B));
 
 %mu=1/2/max(diag(C13yy*D4));
 %B=B+mu*SS;                     % Simbec original.
 
 SSS=SS*B'; SK=(SSS'-SSS)/2;
 mu=mu/2;                        % Simbec stabilized version for noise.
 B=inv(eye(e)+mu*SK/2)*(eye(e)-mu*SK/2)*B+mu*SS;
 
end
  % Convergence?   
    
  phi(1,it+1)=sum(...
     w(2)/(2*alpha(2))*abs(diag(C11yy).^alpha(2))+...
     w(3)/(3*alpha(3))*abs(diag(C12yy).^alpha(3))+...
     w(4)/(4*alpha(4))*abs(diag(C13yy).^alpha(4))+...
     w(5)/(5*alpha(5))*abs(diag(C14yy).^alpha(5))+...
     w(6)/(6*alpha(6))*abs(diag(C15yy).^alpha(6)));   

  delta=(abs(phi(1,it+1)-phi(1,it)))/phi(1,it+1);
  convergence = (delta<stop) | (it==MaxIter);
  
  if (trace(B*B')>1e3)& (reinic>=0), 
     B=rand(e,n);B=B/norm(B);  
     if reinic==0 %| e>1, 
        B=zeros(e,n); convergence=1; 
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
   if verboseKS,
      fprintf('%d  Cum_Index: %6.3f  Convergence: %6.5f\r', it,phi(1,it+1),delta);     
   else
      fprintf('it = %d\n', it );     
 	end     
  
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

C_index=phi(1,it+1);

if verboseKS,
	fprintf('End.');
	fprintf(' \n');
end

% Sort outputs according to their index

[v,i]= sort(w(2)/(2*alpha(2))*abs(diag(C11yy).^alpha(2))+...
            w(3)/(3*alpha(3))*abs(diag(C12yy).^alpha(3))+...
            w(4)/(4*alpha(4))*abs(diag(C13yy).^alpha(4))+...
            w(5)/(5*alpha(5))*abs(diag(C14yy).^alpha(5))+...
            w(6)/(6*alpha(6))*abs(diag(C15yy).^alpha(6)));
        
B=B(flipud(i),:);        