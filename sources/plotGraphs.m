function [] = plotGraphs(s,se,sparse_mix)

[n_sources, ~] = size(s);


for i=1:n_sources
    hFig = figure(i);
    suptitle(strcat('Audio Signal No ',num2str(i)));
    
    subplot(2,2,1);
    h1 = plot(reshape(s(i,:),1,[]),'b');
    legend(h1,{'Original signal'},'FontSize',12);
    axis tight
    
    subplot(2,2,2);
    h2 = plot(se(i,:),'b');
    legend(h2,{'Reconstructed signal'},'FontSize',12);
    axis tight
    
    subplot(2,2,3);
    h2 = plot((s(i,:)-se(i,:)).^2,'r');
    legend(h2,{'MSE'},'FontSize',12);

    axis tight
    set(hFig, 'Position', [0 698 836 170]);
end


%% plotting sparse mixture (1 vs all) 
% for i=1:n_sources-1
%     figure;
%     plot(sparse_mix(1,:),sparse_mix(i+1,:),'.b');
%     axis tight
% end

%% plotting mixes audio sources
% if(view)
%     hFig = figure(num_sources+1);
%     suptitle('Signal Mixture');
% 
%     for i = 1:num_sources
%         subplot(num_sources,1,i);
%         h = plot(X(i,:),'b');
%         legend(h,{strcat('Mixture No ',num2str(i))},'FontSize',10);
%     end
%     set(hFig, 'Position', [0 698 780 170])
% end


end

