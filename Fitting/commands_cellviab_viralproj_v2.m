%% Loading data

[Dead_cell_colour, Tumour_cell_colour, Immune_cell_colour, BiTEs_colour, Virus_colour] = colorscheme();

cell_data = mean([93.34949	111.4028	110.245	    115.3608	112.3317	107.8083	103.0156	100.0135	100.2019;...
            61.16169	68.55299	56.94973	68.89266	72.01087	67.41848	70.26495	66.99049	67.22147;...
            56.29388	63.38899	63.32852	63.54016	65.55582	65.75739	69.00259	74.11227	75.3519;...
            72.01681	70.99303	NaN         69.4616	    70.78151	73.91206	73.56516	78.95479	69.86772]',"omitnan");

cell_data(1) = 100;        
        
cell_time = [24 48 72 96]';

viral_time = [12 24 36 48 72 96];

viral_data = mean([875 	400	    1250;...
12500	9000	2000;...
575000	80000	175000;...
300000	250000	10000;...
2000	4000	350;...
50	    50	    1]');


%parameters
p.r = log(2)/24;%0.004;

p.beta = 1e-2;
p.eta = 100000;
p.d_V = 23/24;
p.alpha = 100;%0.2;%viral_data(1)*0.02;
p.h = 1;

tf = 96;

initialconds = [1e5 1e3];

cell_data_indi = [93.34949	111.4028	110.245	    115.3608	112.3317	107.8083	103.0156	100.0135	100.2019;...
            61.16169	68.55299	56.94973	68.89266	72.01087	67.41848	70.26495	66.99049	67.22147;...
            56.29388	63.38899	63.32852	63.54016	65.55582	65.75739	69.00259	74.11227	75.3519;...
            72.01681	70.99303	NaN         69.4616	    70.78151	73.91206	73.56516	78.95479	69.86772]';

viral_data_indi = [875 	400	    1250;...
12500	9000	2000;...
575000	80000	175000;...
300000	250000	10000;...
2000	4000	350;...
50	    50	    1]';
%%

p.beta = 1;
p.d_V = 0.1;

[time_control,control_U,control_V] = modelsimulator_ODE_case1(p,tf,[initialconds(1),0]);
[time_treat,treat_U,treat_V] = modelsimulator_ODE_case1(p,tf,[initialconds]);

figure
subplot(1,4,1)
hold on 
plot(cell_time,interp1(time_control,control_U,cell_time'))
plot(cell_time, cell_data/100.*interp1(time_control,control_U,cell_time'))
plot(cell_time,interp1(time_treat,treat_U,cell_time'),'--')
legend('Control growth','dataxcontrol','mod sim')
subplot(1,4,2)
plot(time_treat,treat_V)
subplot(1,4,3)
plot(time_treat,treat_V)
set(gca,'yscale','log')
subplot(1,4,4)
plot(cell_time,interp1(time_treat,treat_U,cell_time')./interp1(time_control,control_U,cell_time')*100)



cell_viab = interp1(time_treat,treat_U,cell_time')./interp1(time_control,control_U,cell_time')*100;

Vsol = interp1(time_treat,treat_V,viral_time);
val_cell_viab = [cell_viab- cell_data]./max(cell_data);
val_viral_proj = [Vsol - viral_data./p.alpha]/max(viral_data);
val = [val_cell_viab val_viral_proj];
    sum(abs(val).^2)
        


%% attempting fit of ODE

p.alpha = 30;
p.beta = 1e-6;
p.d_V = 2;

initialconds = [1e5 1E2]; %usualy 1e5

[time, U, V,xmultinonlin,CI] = fit_ODE_model_v2(p,tf,initialconds,cell_data,viral_data,cell_time,viral_time);

%%
[time_control,control_U,control_V] = modelsimulator_ODE_case1(p,tf,[initialconds(1),0]);

figure
hold on 
plot(time,U)
plot(time_control,control_U)
xlabel('Time (hours)')
ylabel('Cell count')
set(gca,'FontSize',18)

figure
hold on 
subplot(1,3,1)
%plot(time, (U)/(initialconds(1))*100,'Color',Tumour_cell_colour,'LineWidth',2)
plot(time, U./control_U*100,'Color',Tumour_cell_colour,'LineWidth',2)
hold on 
plot(cell_time,cell_data,'ko:','LineWidth',2)
xlabel('Time (hours)')
ylabel('Cell viability (%)')

subplot(1,3,2)
plot(time,V,'Color',Virus_colour,'LineWidth',2)
hold on 
plot(viral_time,viral_data,'ko:','LineWidth',2)
set(gca,'yscale','linear')
xlabel('Time (hours)')
ylabel('Viral projeny (linear)')

subplot(1,3,3)
plot(time,V,'Color',Virus_colour,'LineWidth',2)
hold on 
plot(viral_time,viral_data,'ko:','LineWidth',2)
set(gca,'yscale','log')
xlabel('Time (hours)')
ylabel('Viral projeny (log)')

%
cell_data_indi = [93.34949	111.4028	110.245	    115.3608	112.3317	107.8083	103.0156	100.0135	100.2019;...
            61.16169	68.55299	56.94973	68.89266	72.01087	67.41848	70.26495	66.99049	67.22147;...
            56.29388	63.38899	63.32852	63.54016	65.55582	65.75739	69.00259	74.11227	75.3519;...
            72.01681	70.99303	NaN         69.4616	    70.78151	73.91206	73.56516	78.95479	69.86772]';

viral_data_indi = [875 	400	    1250;...
12500	9000	2000;...
575000	80000	175000;...
300000	250000	10000;...
2000	4000	350;...
50	    50	    1]';
%%
figure
hold on 
subplot(1,3,1)
hold on 
plot(cell_time,cell_data_indi,'ko:','LineWidth',1)
plot(time, U/initialconds(1)*100,'Color',Tumour_cell_colour,'LineWidth',2)
xlabel('Time (hours)')
ylabel('Cell viability (%)')
set(gca,'FontSize',18)

subplot(1,3,2)
hold on 
plot(viral_time,viral_data_indi,'ko:','LineWidth',1)
plot(time,V*xmultinonlin(4),'Color',Virus_colour,'LineWidth',2)
set(gca,'yscale','linear')
xlabel('Time (hours)')
ylabel('Viral progeny (linear)')
set(gca,'FontSize',18)

subplot(1,3,3)
hold on 
plot(viral_time,viral_data_indi,'ko:','LineWidth',1)
plot(time,V*xmultinonlin(4),'Color',Virus_colour,'LineWidth',2)
set(gca,'yscale','log')
xlabel('Time (hours)')
ylabel('Viral progeny (log)')
set(gca,'FontSize',18)

% std errorbars
figure
hold on 
subplot(1,3,1)
hold on 
errorbar(cell_time,cell_data,nanstd(cell_data_indi),'k','LineWidth',2)
plot(time, U/initialconds(1)*100,'Color',Tumour_cell_colour,'LineWidth',2)
xlabel('Time (hours)')
ylabel('Cell viability (%)')
set(gca,'FontSize',18)

subplot(1,3,2)
hold on 
errorbar(viral_time,viral_data,nanstd(viral_data_indi),'k','LineWidth',2)
plot(time,V*xmultinonlin(4),'Color',Virus_colour,'LineWidth',2)
set(gca,'yscale','linear')
xlabel('Time (hours)')
ylabel('Viral progeny (linear)')
set(gca,'FontSize',18)

subplot(1,3,3)
hold on 
errorbar(viral_time,viral_data,nanstd(viral_data_indi),'k','LineWidth',2)
plot(time,V*xmultinonlin(4),'Color',Virus_colour,'LineWidth',2)
set(gca,'yscale','log')
xlabel('Time (hours)')
ylabel('Viral progeny (log)')
set(gca,'FontSize',18)

% Simple PL building for beta


Opt = xmultinonlin;

%%

figure
hold on 
subplot(1,2,1)
hold on 
plot(time, U,'Color',Tumour_cell_colour,'LineWidth',2)
xlabel('Time (hours)')
ylabel('Cells (count)')
set(gca,'FontSize',18)

subplot(1,2,2)
hold on 
plot(time,V,'Color',Virus_colour,'LineWidth',2)
set(gca,'yscale','linear')
xlabel('Time (hours)')
ylabel('Virions (count)')
set(gca,'FontSize',18)



% std errorbars
figure
hold on 
errorbar(cell_time,cell_data,nanstd(cell_data_indi),'k','LineWidth',2,'CapSize',12)
plot(time, U/initialconds(1)*100,'Color',Tumour_cell_colour,'LineWidth',2.5)
xlabel('Time (hours)')
ylabel('Cell viability (%)')
set(gca,'FontSize',18)
saveas(gca,'Fig2A.fig')
saveas(gca,'Fig2A.png')

figure
hold on
errorbar(viral_time,viral_data,nanstd(viral_data_indi),'k','LineWidth',2,'CapSize',12)
plot(time,V*xmultinonlin(4),'Color',Virus_colour,'LineWidth',2.5)
set(gca,'yscale','log')
xlabel('Time (hours)')
ylabel('Viral progeny (log)')
set(gca,'FontSize',18)
saveas(gca,'Fig2B.fig')
saveas(gca,'Fig2B.png')
