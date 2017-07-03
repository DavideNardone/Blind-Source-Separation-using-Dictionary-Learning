function [Y,B,varargout]=LUJ1D(X,varargin)
%LU based Jacbi-like JD; This function minimizes the cost
%J_{1}(B)=\sum{i=1}^{N} \|BC_{i}B^{T}-diag(BC_{i}B^{T})\|_{F}^{2}
%where \{C_{i}\}_{i=1}^{N} is a set of N, n\times n symmetric matrices
%and B the joint diagonalizer sought. A related measure that is used
%to measure the error is J_{2}=\sum{i=1}^{N} \|C_{i}-B^{-1}diag(BC_{i}B^{T})B^{-T}\|_{F}^{2}
%
%
%Standard usage: [Y,B]=LUJ1D(X),
%Here X is a large matrix of size n\times nN which contains the
%matrices to be jointly diagonalized such that X=[C1,C2,...,CN],
%Y contains the jointly diagonalized version of the input
%matrices, and B is the found diagonalizer.
%
%
%More controlled usage:[Y,B,S,BB]=LUJ1D(X,'mode',ERR or ITER,RBALANCE):
%
%Inputs:
%'mode'='B' or 'E' or 'N':  In the 'B' mode the stopping criteria at each
%                           step is max(max(abs(LU-I))) which measures
%                           how much the diagonalizer B has changed
%                           after a sweep. In the 'E' mode
%                           the stopping criterion is the difference between
%                           the values of the cost function J2 in two consequtive
%                           updates.In the 'N' mode the stopping criterion is
%                           the number of sweeps over L and U phases.
%
%ERR: In the 'B' mode it specifies the stopping value for the change in B max(max(abs(LU-I))).
%The default value for ERR in this mode and other modes including standard usage
%is ERR=10^-5. In implementation of the algorithm in order to account
%for dpendence of accuracy on the dimension n ERR is multiplied
%by n the size of matrices for JD. In the 'E' mode it ERR specifies the stopping value
%for the relative change of J_{2} in two consequetive sweeps.
%In the 'B' or 'E' mode or the standard mode
%if the change in B or relative change in J2 does not reach ERR after the default number of
%iterations (=200) then the program aborts and itreturns the current computed variables.
%
%ITER: Number of iterations in the 'N' mode
%
%%RBALANCE: if given it is the period for row balancing after each sweep.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%Y= the diagonalized set of matrices
%B=the found joint diagonalizer
%S=a structure containing some information about the run program:
%          S.iterations: number of iterations
%          S.LUerror: the LU error after each sweep
%          S.J2error: the J2 error after each sweep
%          S.J2RelativeError:the relative J2 error after each sweep
%BB=a three dimensional array containing the joint diagonalizer after each sweep
%Note: S and BB are not required as output
%
%This algorithm is based on a paper presented in ICA2006 conference and published in Springer LNCS
%Bijan Afsari, ''Simple LU and QR based Non-Orthogonal Matrix Joint Diagonalization''
%%Coded by Bijan Afsari. Please forward any questions and problem to bijan@glue.umd.edu
%v.1.1
%Acknowledgements: Some data structures and implementation ideas in this code are inspired from the code for JADE
%written by J.F. Cardoso and from the code FFDIAG written by A. Ziehe
%Disclaimer: This code is to be used only for non-commercial research purposes and the author does not
%accept any reponsibility about its performance or fauilure
[n,m]=size(X);N=m/n;

BB=[];
%defaulat values
ERR=1*10^-5;RBALANCE=3;ITER=200;
%%%
MODE='B';
if nargin==0, display('you must enter the data'); B=eye(n); return; end;
if nargin==1, Err=ERR;Rbalance=RBALANCE;end;
if nargin> 1, MODE=upper(varargin{1});
    switch MODE
        case {'B'}
            ERR=varargin{2}; mflag='D'; if ERR >= 1, disp('Error value should be much smaller than unity');B=[];S=[]; return; end;
        case ('E')
            ERR=varargin{2};mflag='E'; if ERR >=1, disp('Error value should be much smaller than unity'); B=[];S=[];return;end;
        case ('N');mflag='N'; ITER=varargin{2}; ERR=0; if ITER <= 1, disp('Number of itternations should be higher than one');B=[];S=[];return;end;
    end
end;
if nargin==4, RBALANCE=varargin{3}; if ceil(RBALANCE)~=RBALANCE | RBALANCE<1, disp('RBALANCE should be a positive integer');B=[];S=[];return;end;end;
JJ=[];EERR=[]; EERRJ2=[];
X1=X;
B=eye(n,n);Binv=eye(n);
J=0;

for t=1:N
    J=J+norm(X1(:,(t-1)*n+1:t*n)-diag(diag(X(:,(t-1)*n+1:t*n))),'fro')^2;
end
JJ=[JJ,J];

%err=10^-3;
%the following part implements a sweep
%%%%%%%%%%%%%%%%%%%%%%%%%%
err=ERR*n+1;
if MODE=='B', ERR=ERR*n;end,
k=0;
while err>ERR & k<ITER
    k=k+1;

    L=eye(n);%Linv=eye(n);
    U=eye(n);%Uinv=eye(n);
    Dinv=eye(n);
    for i=1:n
        for j=i+1:n
            cindex=1:m;
            cindex(i:n:m)=[];
            a=-(X(i,cindex)*X(j,cindex)')/(X(j,cindex)*X(j,cindex)');
            %a=tanh(a);
            if abs(a)>1, a=sign(a)*1;end;
            X(i,:)=a*X(j,:)+X(i,:);
            I=i:n:m;
            J=j:n:m;
            X(:,I)=a*X(:,J)+X(:,I);
            %U=eye(n,n);
            %U(i,j)=a;
            %B=U*B;
            U(i,:)=U(i,:)+a*U(j,:);
            %Uinv(i,:)=Uinv(i,:)-a*Uinv(j,:);
            %Now the L phase

        end%end loop for j
    end

    for i=1:n
        %for j=i+1:n

        rindex=[];
        Xj=[];
        for j=i+1:n
            cindex=1:m;
            cindex(j:n:m)=[];
            a=-(X(i,cindex)*X(j,cindex)')/(X(i,cindex)*X(i,cindex)');
            %coorelation quefficient
            %a=-(X(i,cindex)*X(j,cindex)')/(norm(X(i,cindex))*norm(X(j,cindex)));
            %a=tanh(a);
            if abs(a)>1, a=sign(a)*1; end;
            X(j,:)=a*X(i,:)+X(j,:);
            I=i:n:m;
            J=j:n:m;
            X(:,J)=a*X(:,I)+X(:,J);
            L(j,:)=L(j,:)+a*L(i,:);
            %Linv(j,:)=Linv(j,:)-a*Linv(i,:);
        end%end loop over j
    end


    B=L*U*B;%Binv=Binv*Uinv*Linv;
    %err=norm(L*U-eye(n,n),'fro');
    err=max(max(abs(L*U-eye(n))));EERR=[EERR,err];
    if rem(k,RBALANCE)==0
        d=sum(abs(X'));
        D=diag(1./d*N); Dinv=diag(d*N);
        J=0;
        for t=1:N
            X(:,(t-1)*n+1:t*n)=D*X(:,(t-1)*n+1:t*n)*D;
        end;
        B=D*B; %Binv=Binv*Dinv;
    end
    J=0;
    BB(:,:,k)=B;
    Binv=inv(B);
    for t=1:N
        J=J+norm(X1(:,(t-1)*n+1:t*n)-Binv*diag(diag(X(:,(t-1)*n+1:t*n)))*Binv','fro')^2;
    end
    JJ=[JJ,J];
    if MODE=='E', err=abs(JJ(end-1)-JJ(end))/JJ(end-1);EERRJ2=[EERRJ2,err];end
end
Y=X;
S=struct('iterations',k,'LUerror',EERR,'J2error',JJ,'J2RelativeError',EERRJ2);varargout{1}=S;varargout{2}=BB;