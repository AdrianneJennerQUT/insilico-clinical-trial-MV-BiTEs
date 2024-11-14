function plotting(theta_star,time_grid, U_mat, V_mat, cell_time, cell_data, viral_time, viral_data, distance_vec, B_mat, BITE_time,BITE_data)


figure
subplot(2,2,1)
histogram(theta_star(:,1));
title('\beta')
subplot(2,2,2)
histogram(theta_star(:,2));
title('\eta')
subplot(2,2,3)
histogram(theta_star(:,3));
title('d_V')
subplot(2,2,4)
histogram(theta_star(:,4));
title('\alpha')

figure
plot(distance_vec);

figure
plotmatrix(theta_star);

% simulate all accepted theta_stars


figure
subplot(1,2,1)
hold on 
plot(time_grid,U_mat')
plot(cell_time,cell_data,'ko:','LineWidth',2);
subplot(1,2,2)
hold on
plot(time_grid,V_mat')
plot(viral_time,viral_data,'ko:','LineWidth',2);
set(gca,'yscale','log')


figure 
plot(time_grid,B_mat','LineWidth',2)
hold on 
plot(BITE_time,BITE_data,'ko:','LineWidth',2)
set(gca,'yscale','linear')
xlabel('Time (hours)')
ylabel('Viral projeny (log)')




figure
subplot(5,5,1)
histogram(theta_star(:,1));
subplot(5,5,2)
h = binscatter(theta_star(:,1),theta_star(:,2));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')
subplot(5,5,3)
h = binscatter(theta_star(:,1),theta_star(:,3));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')
subplot(5,5,4)
h = binscatter(theta_star(:,1),theta_star(:,4));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')
subplot(5,5,5)
h = binscatter(theta_star(:,1),theta_star(:,5));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')


subplot(5,5,6)
h = binscatter(theta_star(:,2),theta_star(:,1));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')
subplot(5,5,7)
histogram(theta_star(:,2))
subplot(5,5,8)
h = binscatter(theta_star(:,2),theta_star(:,3));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')
subplot(5,5,9)
h = binscatter(theta_star(:,2),theta_star(:,4));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')
subplot(5,5,10)
h = binscatter(theta_star(:,2),theta_star(:,5));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')


subplot(5,5,11)
h = binscatter(theta_star(:,3),theta_star(:,1));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')
subplot(5,5,12)
h = binscatter(theta_star(:,3),theta_star(:,2));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')
subplot(5,5,13)
histogram(theta_star(:,3))
subplot(5,5,14)
h = binscatter(theta_star(:,3),theta_star(:,4));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')
subplot(5,5,15)
h = binscatter(theta_star(:,3),theta_star(:,5));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')

subplot(5,5,16)
h = binscatter(theta_star(:,4),theta_star(:,1));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')
subplot(5,5,17)
h = binscatter(theta_star(:,4),theta_star(:,2));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')
subplot(5,5,18)
h = binscatter(theta_star(:,4),theta_star(:,4));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')
subplot(5,5,19)
histogram(theta_star(:,4));
subplot(5,5,20)
h = binscatter(theta_star(:,4),theta_star(:,5));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')

subplot(5,5,21)
h = binscatter(theta_star(:,5),theta_star(:,1));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')
subplot(5,5,22)
h = binscatter(theta_star(:,5),theta_star(:,2));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')
subplot(5,5,23)
h = binscatter(theta_star(:,5),theta_star(:,3));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')
subplot(5,5,24)
h = binscatter(theta_star(:,5),theta_star(:,4));
xlim(gca,h.XLimits); 
ylim(gca,h.YLimits); 
h.ShowEmptyBins = 'on';
colormap(gca,'parula')
subplot(5,5,25)
histogram(theta_star(:,5));
end