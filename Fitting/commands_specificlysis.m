%% Loading data
[Dead_cell_colour, Tumour_cell_colour, Immune_cell_colour, BiTEs_colour, Virus_colour] = colorscheme();

BiTE_conc = [1000000 100000 10000 1000 100 0]/1e6;

%Lysis = [24.511943703785153 25.183210093735404 19.85183759861459 7.640396932296106 4.78215453970698 3.5105417960911467];

SD = [21.80327868852459 22.622950819672127 17.04918032786885 4.42622950819672 3.7704918032786843 2.622950819672127];
Lysis = [24.754098360655735 25.24590163934426 20.327868852459012 7.704918032786882 5.409836065573767 4.0983606557377];

%parameters
p.r = 0.004;
p.beta = 11.28;%3.638;
p.eta = 9.0355*1e5;%1.071*1e5;
p.d_V = 0.9628;%1.6017;
p.h = 1;
p.alpha = 100.5547;%135.09;%viral_data(1)*0.02;
p.alpha_2 = 0.001;
p.d_K = 0;
p.gam = 0.01;%1e3; 
p.k= 0.45/24;
p.eps = 0.164;
p.d_B = 0;%-log(0.5)/48;



tf = 48;


%% attempting fit of ODE


[time, U, D, B, K, Lysis_mod xmultinonlin] = fit_ODE_model_Lysis(p,tf,BiTE_conc,Lysis);

time = linspace(0,tf,100);
%%

figure
hold on 
subplot(2,2,1)
plot(time, U,'Color',Tumour_cell_colour,'LineWidth',2)
hold on 
xlabel('Time (hours)')
ylabel('Cells, U(t)')
set(gca,'FontSize',18)

subplot(2,2,2)
plot(time, D,'Color',Dead_cell_colour,'LineWidth',2)
hold on 
xlabel('Time (hours)')
ylabel('Cells, D(t)')
set(gca,'FontSize',18)

subplot(2,2,3)
plot(time,B,'Color',BiTEs_colour,'LineWidth',2)
hold on 
set(gca,'yscale','linear')
xlabel('Time (hours)')
ylabel('BiTEs, B(t)')
set(gca,'FontSize',18)

subplot(2,2,4)
plot(time,K,'Color',Immune_cell_colour,'LineWidth',2)
hold on 
set(gca,'yscale','linear')
xlabel('Time (hours)')
ylabel('Cells, K(t)')
set(gca,'FontSize',18)

figure
hold on 
plot(BiTE_conc,Lysis)
plot(BiTE_conc, Lysis_mod*100,'ko-','LineWidth',2)
xlabel('BiTE')
ylabel('Specific Lysis %')
set(gca,'xscale','log')

figure
hold on 
errorbar([1 2 3 4 5 6],Lysis,Lysis-SD,'k','LineWidth',2,'CapSize',13)
plot([1 2 3 4 5 6], Lysis_mod*100,'Color',Dead_cell_colour,'LineWidth',2.5)
xlabel('BiTE (conc.)')
ylabel('Specific Lysis %')
set(gca,'Xtick',[1 2 3 4 5 6], 'XTickLabelRotation',45,'XTickLabel',{'1\mug/ml','100 ng/ml','10 ng/ml','1 ng/ml','100 pg/ml','0 pg/ml'})
%set(gca,'xscale','log')
set(gca,'FontSize',18)
saveas(gca,'Fig2D.fig')
saveas(gca,'Fig2D.png')



Opt = xmultinonlin;

