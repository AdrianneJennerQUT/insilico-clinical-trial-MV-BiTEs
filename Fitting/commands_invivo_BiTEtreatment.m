
[Dead_cell_colour, Tumour_cell_colour, Immune_cell_colour, BiTEs_colour, Virus_colour] = colorscheme();

Days = [ 5 6 7 13 15 17 20 22 24 27 29 31 34 36];

Tvol = [4.864	31.104	7.0395	29.478	23.328	21.294	12.5	10.368	30.0125	11.9025
        6.6785	27.7695	17.661	58.8245	30.324	18.75	9.702	13.5375	19.22	18.5895
        12.8125	45.63	28.7875	74.0015	37.926	27.048	23.552	30.276	21.4375	41.154
        35	    202	    111	    135	    133	    198	    59	    156	    93	    296
        22	    250	    120	    155	    141	    224	    98	    568	    133	    463
        32	    379	    168	    200	    145	    335	    170	    1598	116	    720
        34	    1263	473	    380	    194		NaN     NaN		NaN     368	    1576
        50		NaN	    NaN     648	    272		NaN     832		NaN     696	    NaN
        64		NaN	    NaN     1180	495		NaN		NaN	    NaN     NaN     NaN
        129		NaN		NaN     NaN     875		NaN		NaN	    NaN     NaN     NaN
        282		NaN 	NaN		NaN		NaN		NaN     NaN     NaN     NaN     NaN
        435		NaN		NaN		NaN		NaN	    NaN     NaN     NaN     NaN     NaN
        712		NaN		NaN		NaN		NaN	    NaN     NaN     NaN     NaN     NaN
        1203	NaN		NaN		NaN		NaN		NaN     NaN     NaN     NaN     NaN];

Mock_col = [94,79,162];
BiTE_col = [50,136,189];
MV_col = [102,194,165];
CEA_col = [171,221,164];
UV_col = [254,224,139];
CD20_col = [244,109,67];

load('invivodata.mat')

sum(isnan(BiTE_vol)');
BiTE_50 = find(sum(isnan(BiTE_vol)')>5);

Tvol_BiTEs = nanmean(BiTE_vol(2:BiTE_50-1,:)');
Days_data = BiTE_time;
Tvol_data = Tvol_BiTEs;

%%
B0 = 2.17*10;
initialconds = [1e5 0 0 B0 2.5*1e5];

% in vivo growth with BiTEs producing OV - assuming immune cells same as in
% vitro experiment

%units days
p.r =0.1321;
p.L = 2.73*1e9;
p.beta = 270.8;
p.h = 1;
p.k = 0.1973;
p.eps = 0.1716;
p.gam = 4.25*1e-3;
p.d_D = 4.8;
p.alpha = 100.56;
p.alpha_2 = 0.001;
p.d_V = 23.11;
p.d_B = 7.89;
p.d_T = 0.35;
p.s = 1;
p.eta = 9.04*1e5;

tf = BiTE_time(end);

% FIT d_B

[time, U, D, V, B, T, xmultinonlin] = fit_ODE_model_BiTEs_invivo(p,tf,initialconds,Days_data,Tvol_data);

%
figure
errorbar(Days(1:7),nanmean(Tvol(1:7,:)'), nanstd(Tvol(1:7,:)'),'k','LineWidth',2,'CapSize',13)
hold on 
plot(time, U/1e6,'Color',Tumour_cell_colour,'LineWidth',2.5)
hold on 
xlabel('Time (days)')
ylabel('Tumour volume (mm^3)')
ylim([0 1200])
box off
set(gca,'FontSize',18)
saveas(gca,'Fig2F.fig')
saveas(gca,'Fig2F.png')


figure
subplot(2,2,1)
plot(time,U)
hold on
plot(time,D)
subplot(2,2,2)
plot(time, B);
legend('B')
subplot(2,2,3)
plot(time, T);
legend('T')
subplot(2,2,4)
plot(time, V);
legend('V')

STOP
%%


figure
p.alpha_2 = 0;
   [time, U, D, V, B, T] = modelsimulator_ODE_fullmod(p,tf,initialconds);
subplot(1,2,1)
hold on 
plot(time, U/1e6,'LineWidth',2)
hold on 
plot(Days_data,Tvol_data,'ko:','LineWidth',2)
xlabel('Time (hours)')
ylabel('Cell viability (%)')

subplot(1,2,2)
plot(time, B);


%%

p.alpha_2 = xmultinonlin(1);
p.s = xmultinonlin(2);

tf = 30;
   [time, U, D, V, B, T] = modelsimulator_ODE_fullmod(p,tf,initialconds);

figure
subplot(1,4,1)
plot(time,U,'Color',Tumour_cell_colour,'LineWidth',3)
hold on
plot(time,D,'Color',Dead_cell_colour,'LineWidth',3)
xlabel('  ')
ylabel('Tumour cells')
legend('U','D')
set(gca,'FontSize',22)
subplot(1,4,2)
plot(time, T,'Color',Immune_cell_colour,'LineWidth',3)
ylabel('CD3 cells (T)')
xlabel('Time (days)')
set(gca,'FontSize',22)
subplot(1,4,3)
plot(time, B,'Color',BiTEs_colour,'LineWidth',3)
xlabel('  ')
ylabel('BiTEs (B)')
set(gca,'FontSize',22)
subplot(1,4,4)
plot(time, V,'Color',Virus_colour,'LineWidth',3)
xlabel('  ')
ylabel('Virus (V)')
set(gca,'FontSize',22)

figure
hold on 
plot(time, U/1e6)
plot(Days,Tvol)

%% simulate varying OV dosage frequency over 60 days


% up to 7 dosages total up to 7 days apart, U(end) measured

initialconds = [1e5 0 1 0 0];%2.5*1e5];
tf = 65;
U_matrix = [];

for totaldosages = 1:9
    for daysapart = 1:7
        
       [U_mat, D_mat, V_mat, B_mat, T_mat, U_end, D_end, V_end, B_end, T_end] = model_sim_varying_sched(p,tf,initialconds,totaldosages,daysapart);

        U_end_matrix(totaldosages,daysapart) = real(U_end);
        D_end_matrix(totaldosages,daysapart) = real(D_end);
        V_end_matrix(totaldosages,daysapart) = real(V_end);
        B_end_matrix(totaldosages,daysapart) = real(B_end);
        T_end_matrix(totaldosages,daysapart) = real(T_end);
        U_FULL_matrix{totaldosages}(daysapart,:) = real(U_mat);
        D_FULL_matrix{totaldosages}(daysapart,:) = real(D_mat);
        V_FULL_matrix{totaldosages}(daysapart,:) = real(V_mat);
        B_FULL_matrix{totaldosages}(daysapart,:)= real(B_mat);
        T_FULL_matrix{totaldosages}(daysapart,:) = real(T_mat);

    end

end

figure
heatmap(real(U_end_matrix))
xlabel('Days apart')
ylabel('Totaldosages')

colormap(jet)


%%
