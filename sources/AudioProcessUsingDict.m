function [w,y] = AudioProcessUsingDict(audio_sample,D,param)
% OUTPUT:
% w is a sparse representation of the audio_sample
% y is the reconstructed representation

x = double(audio_sample);

w = mexOMP(x,D,param);
y = D*w;

end

