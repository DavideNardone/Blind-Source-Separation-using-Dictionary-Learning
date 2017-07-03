close all;
end_res = length(N_SAMPLE);
end_res = end_res-1;

% suptitle('Separation performance in terms of average SDR for four speech sources')
subplot(3,1,1);
h2 = plot(N_SAMPLE(1:end_res), EL_TIME_mean(1:end_res),'b-*','LineWidth',2.5,'markers',16);
xlabel('Block length (samples)','FontSize',30,'Fontweight','bold','Color','black');
ylabel('Time (seconds)','FontSize',30,'Fontweight','bold','Color','black');
set(gca,'Fontsize',33)
ax = gca;
ax.XTick = N_SAMPLE;
axis tight

% suptitle('Required time for separating four speech sources')
subplot(3,1,2)
avg_res = mean(SDR_mean');
i_neg = find(avg_res<0);
plot(N_SAMPLE(1:end_res), avg_res(1:end_res),'b-*','LineWidth',2.5,'markers',16)
hold on
plot(N_SAMPLE(i_neg), avg_res(i_neg),'r-*','LineWidth',2.5,'markers',16)
xlabel('Block length (samples)','FontSize',30,'Fontweight','bold','Color','black');
ylabel('Performance (dB)','FontSize',30,'Fontweight','bold','Color','black');
set(gca,'Fontsize',33) 
ax = gca;
ax.XTick = N_SAMPLE;
axis tight

% suptitle('Average running time/average SDR')
subplot(3,1,3);
h3 = semilogy(N_SAMPLE(1:end_res), EL_TIME_mean(1:end_res)./abs(avg_res(1:end_res)'),'b-*','LineWidth',2.5,'markers',16);
xlabel('Block length (samples)','FontSize',30,'Fontweight','bold','Color','black');
ylabel('Cost-Benefit','FontSize',30,'Fontweight','bold','Color','black');
set(gca,'Fontsize',33)
ax = gca;
ax.XTick = N_SAMPLE;
axis tight


