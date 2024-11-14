[Dead_cell_colour, Tumour_cell_colour, Immune_cell_colour, BiTEs_colour, Virus_colour] = colorscheme()


cell_data = nanmean([93.34949	111.4028	110.245	    115.3608	112.3317	107.8083	103.0156	100.0135	100.2019;...
            61.16169	68.55299	56.94973	68.89266	72.01087	67.41848	70.26495	66.99049	67.22147;...
            56.29388	63.38899	63.32852	63.54016	65.55582	65.75739	69.00259	74.11227	75.3519;...
            72.01681	70.99303	NaN         69.4616	    70.78151	73.91206	73.56516	78.95479	69.86772]');

cell_time = [24 48 72 96]';

viral_time = [12 24 36 48 72 96];

viral_data = mean([875 	400	    1250;...
12500	9000	2000;...
575000	80000	175000;...
300000	250000	10000;...
2000	4000	350;...
50	    50	    1]');

%% Loading data

%BITE_data = [3.177218764912695 54.78481900254186 59.09954129780583 59.64803615079714]*10;
BITE_data = [1.468115942 2.424637681 4.155072464 1.791304348];
BITE_time = [24 48 72 96];


%parameters
p.r = 0.004;
p.beta = 11.28;%3.638;
p.eta = 9.0355*1e5;%1.071*1e5;
p.d_V = 0.9628;%1.6017;
p.h = 1;
p.alpha = 100.5547;%135.09;%viral_data(1)*0.02;
p.alpha_2 = 2;%1e-3;


tf = 96;


%% attempting fit of ODE

initialconds = [1e5 1 0];

[time, U, V, B, xmultinonlin] = fit_ODE_model_BiTEs(p,tf,initialconds,BITE_data,BITE_time);


%
figure
hold on 
subplot(1,3,1)
plot(time, (U)/(initialconds(1))*100,'Color',Tumour_cell_colour,'LineWidth',2)
hold on 
plot(cell_time,cell_data,'ko:','LineWidth',2)
xlabel('Time (hours)')
ylabel('Cell viability (%)')
set(gca,'FontSize',18)
box off

subplot(1,3,2)
plot(time,V*p.alpha,'Color',Virus_colour,'LineWidth',2)
hold on 
plot(viral_time,viral_data,'ko:','LineWidth',2)
set(gca,'yscale','log')
xlabel('Time (hours)')
ylabel('Viral projeny (linear)')
set(gca,'FontSize',18)
box off

subplot(1,3,3)
plot(time,B,'Color',BiTEs_colour,'LineWidth',2)
hold on 
plot(BITE_time,BITE_data,'ko:','LineWidth',2)
set(gca,'yscale','linear')
xlabel('Time (hours)')
ylabel('Viral projeny (log)')
set(gca,'FontSize',18)
box off

% Simple PL building for beta

figure
plot(time,B,'Color',BiTEs_colour,'LineWidth',2.5)
hold on 
plot(BITE_time,BITE_data,'ko-','LineWidth',2)
set(gca,'yscale','linear')
xlabel('Time (hours)')
ylabel('BiTEs (fold change)')
set(gca,'FontSize',18)
box off
saveas(gca,'Fig2C.fig')
saveas(gca,'Fig2C.png')


Opt = xmultinonlin;

%%

