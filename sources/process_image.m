function [D,w,y2]=process_image(image_data)

% image_data=imread(filename); %Read data from image
orig=double(image_data);
no_of_images=size(image_data(1,:),2); %No. of columns in x(equivalent to number of signals)
dim=size(image_data(:,1),1); %No. of rows in x (equivalent to dimension of each signal)
subdiv=8; % Size of subdivision to compute sparse representation
nob=(dim/subdiv)*(no_of_images/subdiv); % Total number of 8x8 Blocks

%% Dividing the image into 8x8 Blocks (patches)
Block=[];
kk=0;
for i=1:(dim/subdiv)
    for j=1:(no_of_images/subdiv)
        Block(:,:,kk+j)=orig((subdiv*(i-1)+1:subdiv*(i-1)+subdiv),(subdiv*(j-1)+1:subdiv*(j-1)+subdiv));
    end
    kk=kk+(dim/subdiv);
end

%% Dictionary Learning
d_size=256; %Size of dictionary (No. of atoms)
no_dict=floor(d_size/subdiv);
r1=randi(nob,no_dict,1);
D1=Block(:,:,r1);
A=D1(:,:,1);
for i=2:size(D1,3)
    A=horzcat(A,D1(:,:,i));
end
D=A;
D=mexNormalize(D);
x=reshape(orig,subdiv,[]);
noIt=50; %No. of iterations for learning
param.L=2; %Maximum 0-norm of each column in sparse representation
  
% Method of Optimal Directions (MOD)
for it=1:noIt
    w=mexOMP(x,D,param); %Finds a sparse representaion of w using mexOMP
    xwt=mexCalcXAt(x,w);
    wwt=mexCalcAAt(w);
    ww_inv=mexInvSym(wwt);
    D=mexCalcXY(xwt,ww_inv);
    D=mexNormalize(D); %Make each column have unit norm
end
w=mexOMP(x,D,param);
y=D*w;
y2=reshape(y,dim,no_of_images);


end