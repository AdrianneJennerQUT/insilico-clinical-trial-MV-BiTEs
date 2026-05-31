%% initialise parameters

%load colour scheme
[Dead_cell_colour, Tumour_cell_colour, Immune_cell_colour, BiTEs_colour, Virus_colour] = colorscheme();

%load parameters
p = parameters_invivo();

tf = 25;

initialconds = [p.U0 0 p.V0 p.B0 p.K0];

[time,U,D,V,B,T] = fullmodel(p,initialconds,tf);

%%
figure
subplot(1,4,1)
hold on 
plot(time,U,'Color',Tumour_cell_colour,'LineWidth',2)
plot(time,D,'Color',Dead_cell_colour,'LineWidth',2)
xlabel('Time (days)')
ylabel('Tumour cells')
legend('U','D')
set(gca,'FontSize',18)


subplot(1,4,2)
hold on 
plot(time,V,'Color',Virus_colour,'LineWidth',2)
xlabel('Time (days)')
ylabel('Virus (V)')
set(gca,'FontSize',18)


subplot(1,4,3)
hold on 
plot(time,B,'Color',BiTEs_colour,'LineWidth',2)
xlabel('Time (days)')
ylabel('BiTEs (B)')
set(gca,'FontSize',18)


subplot(1,4,4)
hold on 
plot(time,T,'Color',Immune_cell_colour,'LineWidth',2)
xlabel('Time (days)')
ylabel('CD3 cells (T)')
set(gca,'FontSize',18)

