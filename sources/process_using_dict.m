function [w,y2]=process_using_dict(image_data,D)

orig=double(image_data);
no_of_images=size(image_data(1,:),2); %No. of columns in x(equivalent to number of signals)
dim=size(image_data(:,1),1);          %No. of rows in x (equivalent to dimension of each signal)
subdiv=8;                       % Size of subdivision to compute sparse representation
x=reshape(orig,subdiv,[]);
%     d_size=256;          %Size of dictionary(No. of atoms)
%     r=randi(no_of_images,d_size,1);
%     D=x(:,r);          %Initialize Dictionary with random columns from x
%    D=mexNormalize(D); %Make each column have unit norm
%noIt=5;   %No. of iterations for learning
param.L=2; %Maximum 0-norm of each column in sparse representation

%for it=1:noIt
w=mexOMP(x,D,param);    %Finds a sparse representaion of w using mexOMP
%             xwt=mexCalcXAt(x,w);
%             wwt=mexCalcAAt(w);
%             ww_inv=mexInvSym(wwt);
%             D=mexCalcXY(xwt,ww_inv);
%             D=mexNormalize(D);  %Make each column have unit norm
%end
y=D*w;
y2=reshape(y,dim,no_of_images);

