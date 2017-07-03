clear all
close all
clc
% run C:\Users\Shreshth\Downloads\spams-matlab-v2.5-svn2014-07-04\spams-matlab\start_spams.m
subdiv=8;                       % Size of subdivision to compute sparse representation
d_size=256;          %Size of dictionary(No. of atoms)

[audio_data1,fs1]=audioread('/Users/davidenardone/Dropbox/MATLAB/SM/dataset/AUDIO/cmu_us_ked_timit/wav/kdt_335.wav');
[D1, w1, y1]=process_audio(audio_data1);

[audio_data2,fs2]=audioread('/Users/davidenardone/Dropbox/MATLAB/SM/dataset/AUDIO/cmu_us_ked_timit/wav/kdt_349.wav');
[w2, y2]=process_using_dict(audio_data2,D1);

w1=reshape(w1,1,[]);
w2=reshape(w2,1,[]);
x1=reshape(double(audio_data1),1,[]);
x2=reshape(double(audio_data2),1,[]);

A=abs(randn(2,2));
A=A./repmat(sqrt(sum(A.^2)),[size(A,1) 1]);

x=[x1;x2];
y_mix=A*x;
y_final1=reshape(y_mix(1,:),1,[]);

audiowrite('mixes/mix1.wav',y_mix(1,:),fs1);
y_final2=reshape(y_mix(2,:),1,[]);

audiowrite('mixes/mix2.wav',y_mix(2,:),fs2);


[wm1,ym1]=process_using_dict(y_final1,D1);
wm1=reshape(wm1,1,[]);
[wm2,ym1]=process_using_dict(y_final2,D1);
wm2=reshape(wm2,1,[]);
wm=[wm1;wm2];

param.L=2;
w_rec=zeros(2,size(wm,2));
[i,j,v]=find(wm);
for i=1:numel(j)
    w_rec(:,j(i))=mexOMP(full(wm(:,j(i))),A,param);
end
s1=reshape(w_rec(1,:),256,[]);
s1=reshape(D1*s1,1,[]);
audiowrite('output/extract1.wav',s1,fs1);

s2=reshape(w_rec(2,:),256,[]);
s2=reshape(D1*s2,1,[]);
audiowrite('output/extract2.wav',s2,fs2);