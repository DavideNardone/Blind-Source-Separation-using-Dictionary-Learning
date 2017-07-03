function Y = splitaudio(F,len_frame)


%% splitting audio in frames
%         len_frame = 256;
n_frame = ceil(length(F)/len_frame);


Y = zeros(len_frame,n_frame);

for k=1:size(Y,2)
    
    
    startIdx = (len_frame*(k-1))+1;
    endIdx = (len_frame*k);
    
    if(endIdx > length(F))
        F=padarray(F,abs(length(F)-endIdx),'post');
    end
    %             [startIdx endIdx]
    Y(:,k) = F(startIdx:1:endIdx);
end



end

