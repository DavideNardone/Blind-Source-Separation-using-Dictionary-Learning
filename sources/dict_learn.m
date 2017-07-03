clear all
close all
clc
% run C:\Users\Shreshth\Downloads\spams-matlab-v2.5-svn2014-07-04\spams-matlab\start_spams.m
subdiv=8;                       % Size of subdivision to compute sparse representation
d_size=256;          %Size of dictionary(No. of atoms)
image_data1=imread('lena.png'); %Read data from image
[D1, w1, y1]=process_image(image_data1); %Find sparse representation and learn dictionary D1 from image itself

image_data2=imread('zelda.bmp'); %Read data from image
[w2, y2]=process_using_dict(image_data2,D1);    %Find sparse representation of 2nd image in D1

w1=reshape(w1,1,[]);
w2=reshape(w2,1,[]);


x1=reshape(double(image_data1),1,[]);   %Reshaping image as vector for mixing
x2=reshape(double(image_data2),1,[]);   
A=abs(randn(2,2));                      %Random Mixing Vector A
A=A./repmat(sqrt(sum(A.^2)),[size(A,1) 1]); %Normalizing Columns of A
x=[x1;x2];
y_mix=A*x;                              %Forming the mixtures


y_final1=reshape(y_mix(1,:),512,512);
imwrite(uint8(y_final1),'/Users/davidenardone/Dropbox/MATLAB/SM/spams-matlab/data/mix/mix1.bmp','bmp');

y_final2=reshape(y_mix(2,:),512,512);
imwrite(uint8(y_final2),'/Users/davidenardone/Dropbox/MATLAB/SM/spams-matlab/data/mix/mix2.bmp','bmp');


[wm1,ym1]=process_using_dict(y_final1,D1);  %Sparsifying the mixtures
wm1=reshape(wm1,1,[]);
[wm2,ym1]=process_using_dict(y_final2,D1);
wm2=reshape(wm2,1,[]);
wm=[wm1;wm2];
[i,j,v]=find(wm);
wm2=wm(:,j);
[idx,C]=kmeans(full(wm2'),2);                     %Estimating the mixing matrix
C=full(C);
C=C./repmat(sqrt(sum(C.^2)),[size(C,1) 1]);

param.L=2;
w_rec=zeros(2,size(wm,2));
% [i,j,v]=find(wm);
for i=1:numel(j)
    w_rec(:,j(i))=mexOMP(full(wm(:,j(i))),C,param); %Separating the sources
end
s1=reshape(w_rec(1,:),256,[]);
s1=reshape(D1*s1,512,512);
imwrite(uint8(s1),'/Users/davidenardone/Dropbox/MATLAB/SM/spams-matlab/data/extract/extract1.bmp','bmp');

s2=reshape(w_rec(2,:),256,[]);
s2=reshape(D1*s2,512,512);
imwrite(uint8(s2),'/Users/davidenardone/Dropbox/MATLAB/SM/spams-matlab/data/extract/extract2.bmp','bmp');

    