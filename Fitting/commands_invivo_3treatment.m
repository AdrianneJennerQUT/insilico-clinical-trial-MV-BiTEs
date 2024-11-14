%% load data
[Dead_cell_colour, Tumour_cell_colour, Immune_cell_colour, BiTEs_colour, Virus_colour] = colorscheme();

Mock_col = [94,79,162];
BiTE_col = [50,136,189];
MV_col = [102,194,165];
CEA_col = [171,221,164];
UV_col = [254,224,139];
CD20_col = [244,109,67];

load('invivodata.mat')

CD20_vol_NEW = CD20_vol(:,[1:7,10]);

sum(isnan(Mock_vol)');
Mock_50 = find(sum(isnan(Mock_vol)')>5);

sum(isnan(BiTE_vol)');
BiTE_50 = find(sum(isnan(BiTE_vol)')>5);

sum(isnan(MV_vol)');
MV_50 = find(sum(isnan(MV_vol)')>5);

sum(isnan(CD20_vol)');
CD20_50 = find(sum(isnan(CD20_vol)')>5);

Tvol_BiTEs = nanmean(BiTE_vol(2:BiTE_50-1,:)');
Tvol_MV = nanmean(MV_vol(2:MV_50-1,:)');
Tvol_CD20 = nanmean(CD20_vol_NEW(2:CD20_50-1,:)'); %CHANGED TO NO M8 AND M9



%% load parameters

%units days
p.r =0.1321;
p.L = 2.73*1e9;
p.beta = 0.8;%270.8; % XXXXXXXXXXXXXX
p.h = 1;
p.k = 0.1973;
p.eps = 0.1716;
p.gamma = 4.25*1e-3;
p.d_D = 4.8; % XXXXXXXXXXXXXX
p.alpha = 100.56;
p.alpha_2 = 150;%0.001;
p.d_V = 23.11;
p.d_B = 7.89; % XXXXXXXXXXXXXX
p.d_T = 0.35;
p.s = 5; % XXXXXXXXXXXXXX
p.eta = 9.04*1e5;

p.B0 = 21.7;%[1E3]; 
p.V0 = 1e6;%[1e3]; 
p.K0 = 2.5*1e5;%1e3;
p.U0 = 1e6;

tf = 30;

    p.s = 7.6786;%param(1);
    p.d_D = 0.4584;%param(2);
    p.d_B = 0.3825;%param(5);
    p.B0 = 1e6; % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  
%% TESTING MODEL

  
    initialconds = [1e5 0 p.B0 0];%p.K0];
    [timeB, modB] = modelsimulator_ODE_case5(p,tf,initialconds); %BiTEs
    initialconds = [1e5 0 p.V0 p.K0];
    [timeMV, modMV] = modelsimulator_ODE_case6(p,tf,initialconds); %MV
    initialconds = [1e5 0 p.V0 p.B0 p.K0];
    [timeMVB, modMVB] = modelsimulator_ODE_case7(p,tf,initialconds); %MV-BiTEs
    

%%
figure
errorbar(BiTE_time(2:BiTE_50-1),Tvol_BiTEs,nanstd(BiTE_vol(2:BiTE_50-1,:)'),'k','LineWidth',2,'CapSize',13)
hold on 
plot(timeB, modB(:,1)/1e6,'Color',Tumour_cell_colour,'LineWidth',2.5)
hold on 
xlabel('Time (days)')
ylabel('Tumour volume (mm^3)')
ylim([0 1200])
box off
set(gca,'FontSize',18)
saveas(gca,'Fig2F.fig')
saveas(gca,'Fig2F.png')

%
figure
errorbar(MV_time(2:MV_50-1),Tvol_MV,nanstd(MV_vol(2:MV_50-1,:)'),'k','LineWidth',2,'CapSize',13)
hold on 
plot(timeMV, modMV(:,1)/1e6,'Color',Tumour_cell_colour,'LineWidth',2.5)
hold on 
xlabel('Time (days)')
ylabel('Tumour volume (mm^3)')
ylim([0 1200])
box off
set(gca,'FontSize',18)
saveas(gca,'Fig2F.fig')
saveas(gca,'Fig2F.png')


%
figure
errorbar(CD20_time(2:CD20_50-1),Tvol_CD20,nanstd(CD20_vol(2:CD20_50-1,:)'),'k','LineWidth',2,'CapSize',13)
hold on 
plot(timeMVB, modMVB(:,1)/1e6,'Color',Tumour_cell_colour,'LineWidth',2.5)
hold on 
xlabel('Time (days)')
ylabel('Tumour volume (mm^3)')
ylim([0 1200])
box off
set(gca,'FontSize',18)
saveas(gca,'Fig2F.fig')
saveas(gca,'Fig2F.png')

figure
hold on 
errorbar(Mock_time(1:Mock_50-1),nanmean(Mock_vol(1:Mock_50-1,:)'),nanstd(Mock_vol(1:Mock_50-1,:)'),'.-','Color',Mock_col/255,'MarkerSize',20,'LineWidth',2)
errorbar(BiTE_time(2:BiTE_50-1),Tvol_BiTEs,nanstd(BiTE_vol(2:BiTE_50-1,:)'),'.-','Color',BiTE_col/255,'MarkerSize',20,'LineWidth',2)
errorbar(MV_time(2:MV_50-1),Tvol_MV,nanstd(MV_vol(2:MV_50-1,:)'),'.-','Color',MV_col/255,'MarkerSize',20,'LineWidth',2)
errorbar(CD20_time(2:CD20_50-1),Tvol_CD20,nanstd(CD20_vol(2:CD20_50-1,:)'),'.-','Color',CD20_col/255,'MarkerSize',20,'LineWidth',2)

plot(timeB, modB(:,1)/1e6,'k','LineWidth',1.5)
plot(timeMV, modMV(:,1)/1e6,'k','LineWidth',1.5)
plot(timeMVB, modMVB(:,1)/1e6,'k','LineWidth',1.5)
xlim([0 30])
ylim([0 1200])
xlabel('Time (days)')
ylabel('Tumour volume (mm^3)')
ylim([0 1200])
box off
set(gca,'FontSize',18)
legend('Mock','BiTEs','MV','MV-BiTEs')


%% fit model, fit s, alpha_B, beta, d_D
[timeB, modB,timeMV, modMV, timeMVB, modMVB, xmultinonlin] = fit_ODE_model_3treat_invivo(p,tf,...
    Tvol_BiTEs,Tvol_MV,Tvol_CD20,BiTE_time(2:BiTE_50-1),MV_time(2:MV_50-1),CD20_time(2:CD20_50-1));

%%
figure
errorbar(BiTE_time(2:BiTE_50-1),Tvol_BiTEs,nanstd(BiTE_vol(2:BiTE_50-1,:)'),'k','LineWidth',2,'CapSize',13)
hold on 
plot(timeB, modB(:,1)/1e6,'Color',Tumour_cell_colour,'LineWidth',2.5)
hold on 
xlabel('Time (days)')
ylabel('Tumour volume (mm^3)')
ylim([0 1200])
xlim([0 max(BiTE_time(2:BiTE_50-1))])
box off
set(gca,'FontSize',18)
saveas(gca,'Fig2F.fig')
saveas(gca,'Fig2F.png')

%
figure
errorbar(MV_time(2:MV_50-1),Tvol_MV,nanstd(MV_vol(2:MV_50-1,:)'),'k','LineWidth',2,'CapSize',13)
hold on 
plot(timeMV, modMV(:,1)/1e6,'Color',Tumour_cell_colour,'LineWidth',2.5)
hold on 
xlabel('Time (days)')
ylabel('Tumour volume (mm^3)')
ylim([0 1200])
xlim([0 max(MV_time(2:MV_50-1))])
box off
set(gca,'FontSize',18)
saveas(gca,'Fig2F.fig')
saveas(gca,'Fig2F.png')


%
figure
errorbar(CD20_time(2:CD20_50-1),Tvol_CD20,nanstd(CD20_vol(2:CD20_50-1,:)'),'k','LineWidth',2,'CapSize',13)
hold on 
plot(timeMVB, modMVB(:,1)/1e6,'Color',Tumour_cell_colour,'LineWidth',2.5)
hold on 
xlabel('Time (days)')
ylabel('Tumour volume (mm^3)')
ylim([0 1200])
xlim([0 max(CD20_time(2:CD20_50-1))])
box off
set(gca,'FontSize',18)
saveas(gca,'Fig2F.fig')
saveas(gca,'Fig2F.png')



%%
figure
subplot(3,5,1)
plot(timeB,modB(:,1))
legend('U')
subplot(3,5,2)
plot(timeB,modB(:,2));
legend('D')
subplot(3,5,4)
plot(timeB,modB(:,3));
legend('B')
subplot(3,5,5)
plot(timeB, modB(:,4));
legend('T')

subplot(3,5,6)
plot(timeMV,modMV(:,1))
legend('U')
subplot(3,5,7)
plot(timeMV,modMV(:,2));
legend('D')
subplot(3,5,8)
plot(timeMV,modMV(:,3));
legend('V')
subplot(3,5,10)
plot(timeMV, modMV(:,4));
legend('T')


subplot(3,5,11)
plot(timeMVB,modMVB(:,1))
legend('U')
subplot(3,5,12)
plot(timeMVB,modMVB(:,2));
legend('D')
subplot(3,5,13)
plot(timeMVB,modMVB(:,3));
legend('V')
subplot(3,5,14)
plot(timeMVB,modMVB(:,4));
legend('B')
subplot(3,5,15)
plot(timeMVB, modMVB(:,5));
legend('T')

%%
figure
subplot(1,4,1)
hold on 
plot(timeMVB,modMVB(:,1),'Color',Tumour_cell_colour,'LineWidth',2)
plot(timeMVB,modMVB(:,2),'Color',Dead_cell_colour,'LineWidth',2)
xlabel('Time (days)')
ylabel('Tumour cells')
legend('U','D')
set(gca,'FontSize',18)


subplot(1,4,2)
hold on 
plot(timeMVB,modMVB(:,3),'Color',Virus_colour,'LineWidth',2)
xlabel('Time (days)')
ylabel('Virus (V)')
set(gca,'FontSize',18)


subplot(1,4,3)
hold on 
plot(timeMVB,modMVB(:,4),'Color',BiTEs_colour,'LineWidth',2)
xlabel('Time (days)')
ylabel('BiTEs (B)')
set(gca,'FontSize',18)


subplot(1,4,4)
hold on 
plot(timeMVB,modMVB(:,5),'Color',Immune_cell_colour,'LineWidth',2)
xlabel('Time (days)')
ylabel('CD3 cells (T)')
set(gca,'FontSize',18)








