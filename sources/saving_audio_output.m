function [] = saving_audio_output(mix,se,fs,mixPath,outputPath)


[num_sources,~] = size(se);

%% Saving mixtures

% Create result folder for reconstructed sources (if it doesnt exits)
if exist(mixPath,'dir') > 0
    rmdir(mixPath,'s'); %doing so i make sure to have deleted all file within the folder
end
mkdir(mixPath);

for i = 1:num_sources
    mixSourceName = strcat(mixPath,'mix',num2str(i),'.wav');
    audiowrite(mixSourceName,mix(i,:),fs);
end

%% Saving reconstruct audio signals

% Create result folder for reconstructed sources (if it doesnt exits)
if exist(outputPath,'dir') > 0
    rmdir(outputPath,'s'); %doing so i make sure to have deleted all file within the folder
end
mkdir(outputPath);

%saving output
for i = 1:num_sources
    outExPathName = strcat(outputPath,'extract',num2str(i),'.wav');
    audiowrite(outExPathName,se(i,:),fs);
end

end

