%% Loading data
[Dead_cell_colour, Tumour_cell_colour, Immune_cell_colour, BiTEs_colour, Virus_colour] = colorscheme()

Days = [5 6 7 13 15 17];

Tvol = [4.335	26.4275	5.2345	10.625	17.496	6.2135	17	    12.1945	21.97	11.638
        6.615	44.712	12.506	11.232	26.656	10.584	22.646	20.412	25.35	25.047
        22.528	95.616	8.1585	41.8275	30.184	40.3855	17.8605	29.403	26.4915	35.2
        259		NaN     120	    402	    331	    945	    218	    428	    673	    807
        NaN		NaN     223		NaN     571		NaN     474	    310	    1053	1212
        NaN	    NaN     750		NaN     915		NaN     1183	891		NaN     NaN];

Tvol_mean = nanmean(Tvol');

figure
hold on 
plot(Days,Tvol_mean,'k*-')
plot(Days,Tvol,'o--')

tf = 17;

p.L = 1e10;
p.r = 1;


%% attempting fit of ODE


[time, U xmultinonlin] = fit_ODE_model_invivo(p,tf,Days,Tvol_mean);


figure
hold on 
subplot(1,2,1)
plot(time, U/1e6,'Color',Tumour_cell_colour,'LineWidth',2.5)
hold on 
plot(Days,Tvol_mean,'k*-')
xlabel('Time (hours)')
ylabel('Cells, U(t)')

subplot(1,2,2)
plot(time, U/1e6,'Color',Tumour_cell_colour,'LineWidth',2.5)
hold on 
plot(Days,Tvol,'o--')
hold on 
xlabel('Time (hours)')
ylabel('Cells, D(t)')

figure
hold on 
errorbar(Days, nanmean(Tvol'),nanstd(Tvol'),'k','LineWidth',2, 'CapSize',13)
plot(time, U/1e6,'Color',Tumour_cell_colour,'LineWidth',2.5)
ylabel('Tumour volume (mm^3)')
xlabel('Time (days)')
set(gca,'FontSize',18)
xlim([0 18])
set(gca,'XTick',[0 6, 12, 18])
saveas(gca,'Fig2E.fig')
saveas(gca,'Fig2E.png')


Opt = xmultinonlin;

