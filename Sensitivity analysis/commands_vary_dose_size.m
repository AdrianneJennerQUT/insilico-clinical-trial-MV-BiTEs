%% determine impact of varying V0 and B0 size for fixed dosage protocol


%load parameters
p = parameters_invivo();

tf = 25;


%%

dose_V = 1;
dose_B = 1;
        
initial_conditions = [p.U0 0 p.V0*dose_V p.B0*dose_B p.K0];

[time,U,D,V,B,T] = fullmodsim_dosagevary(p,initial_conditions,tf,dose_V,dose_B);
        
time_s = time(2:end);
baseline_AUC  = sum(U(1:end-1).*[time_s-time(1:end-1)]);

dosage_vary = linspace(1e-5,10,100);

for i = 1:length(dosage_vary)
    
    for j = 1:length(dosage_vary)
        
        dose_V = dosage_vary(i);
        dose_B = dosage_vary(j);
        
        initial_conditions = [p.U0 0 p.V0*dose_V p.B0*dose_B p.K0];

        [time,U,D,V,B,T] = fullmodsim_dosagevary(p,initial_conditions,tf,dose_V,dose_B);
        
        time_s = time(2:end);
        AUC_U(i,j) = sum(U(1:end-1).*[time_s-time(1:end-1)])-baseline_AUC;
        U_vec(i,j) = U(end);

    end
    i
    
end
%%

figure
h = heatmap(real(AUC_U),'GridVisible','off')
colm = cbrewer2('PuBu');
colormap(colm);
xlabel('B_0 multipler')
ylabel('V_0 multipler')
set(gca,'FontSize',18)

XLabels = 1:100;
% Convert each number in the array into a string
CustomXLabels = string(XLabels);
% Replace all but the fifth elements by spaces
CustomXLabels(mod(XLabels,100) ~= 0) = " ";

CustomXLabels(1) = '10^{-5}';
CustomXLabels(10) = '1';
CustomXLabels(20) = '2';
CustomXLabels(30) = '3';
CustomXLabels(40) = '4';
CustomXLabels(50) = '5';
CustomXLabels(60) = '6';
CustomXLabels(70) = '7';
CustomXLabels(80) = '8';
CustomXLabels(90) = '9';
CustomXLabels(100) = '10';

h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = CustomXLabels;

%%


dose_V = 1e-5;
dose_B = 1;
        
initial_conditions = [p.U0 0 p.V0*dose_V p.B0*dose_B p.K0];

[time_1eneg5,U_1eneg5,D,V_1eneg5,B,T_1eneg5] = fullmodsim_dosagevary(p,initial_conditions,tf,dose_V,dose_B);



dose_V = 1e-2;
dose_B = 1;
        
initial_conditions = [p.U0 0 p.V0*dose_V p.B0*dose_B p.K0];

[time_1eneg2,U_1eneg2,D,V_1eneg2,B,T_1eneg2] = fullmodsim_dosagevary(p,initial_conditions,tf,dose_V,dose_B);




dose_V = 1;
dose_B = 1;
        
initial_conditions = [p.U0 0 p.V0*dose_V p.B0*dose_B p.K0];

[time_1,U_1,D,V_1,B,T_1] = fullmodsim_dosagevary(p,initial_conditions,tf,dose_V,dose_B);



dose_V = 2;
dose_B = 1;
        
initial_conditions = [p.U0 0 p.V0*dose_V p.B0*dose_B p.K0];

[time_2,U_2,D,V_2,B,T_2] = fullmodsim_dosagevary(p,initial_conditions,tf,dose_V,dose_B);

figure
hold on 
plot(time_1eneg5,U_1eneg5,'Color',[0,104,55]/255,'LineWidth',2.5)
plot(time_1eneg2,U_1eneg2,'Color',[49,163,84]/255,'Linewidth',2.5)
plot(time_1,U_1,'Color',[120,198,121]/255,'LineWidth',2.5)
plot(time_2,U_2,'Color',[194,230,153]/255,'LineWidth',2.5)
legend('10^{-5}','10^{-2}','1','2')
xlabel('Time (days)')
ylabel('Tumour cells U(t)')
set(gca,'FontSize',18)

figure
hold on 
plot(time_2,V_2,'Color',[194,230,153]/255,'LineWidth',2.5)
plot(time_1,V_1,'Color',[120,198,121]/255,'LineWidth',2.5)
plot(time_1eneg2,V_1eneg2,'Color',[49,163,84]/255,'Linewidth',2.5)
plot(time_1eneg5,V_1eneg5,'Color',[0,104,55]/255,'LineWidth',2.5)
legend('2','1','10^{-2}','10^{-5}')
xlabel('Time (days)')
ylabel('Virus V(t)')
set(gca,'FontSize',18)

figure
hold on 
plot(time_1eneg5,T_1eneg5,'Color',[0,104,55]/255,'LineWidth',2.5)
plot(time_1eneg2,T_1eneg2,'Color',[49,163,84]/255,'Linewidth',2.5)
plot(time_1,T_1,'Color',[120,198,121]/255,'LineWidth',2.5)
plot(time_2,T_2,'Color',[194,230,153]/255,'LineWidth',2.5)
legend('10^{-5}','10^{-2}','1','2')
xlabel('Time (days)')
ylabel('T cells T(t)')
set(gca,'FontSize',18)



