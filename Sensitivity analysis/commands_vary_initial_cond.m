% varying initial conditions



%load parameters
p = parameters_invivo();

tf = 25;

%%

cancer_array = [1e1 1e2 1e3 1e4 1e5 1e6 1e7];
immune_array = cancer_array;

for initial_cancer = 1:length(cancer_array)
    
    for initial_immune = 1:length(immune_array)
        
        p.U0 = cancer_array(initial_cancer);
        p.K0 = immune_array(initial_immune);
        
        initial_conditions = [p.U0 0 p.V0 p.B0 p.K0];

        [time,U,D,V,B,T] = fullmodsim(p,initial_conditions,tf);
    
        time_s = time(2:end);
        AUC_U(initial_cancer,initial_immune) = sum(U(1:end-1).*[time_s-time(1:end-1)]);
        U_vec(initial_cancer,initial_immune) = U(end);
    
    end
    initial_cancer
end



figure
h = heatmap(real(U_vec),'GridVisible','off','CellLabelColor','none')
colm = cbrewer2('PuBu');
colormap(colm);
xlabel('T_0 (log)')
ylabel('U_0 (log)')
set(gca,'FontSize',18)

%%

cancer_array = [1e1 1e2 1e3 1e4 1e5 1e6 1e7];
virus_array = cancer_array;

for initial_cancer = 1:length(cancer_array)
    
    for initial_virus = 1:length(virus_array)
        
        p.U0 = cancer_array(initial_cancer);
        p.V0 = immune_array(initial_virus);
        
        initial_conditions = [p.U0 0 p.V0 p.B0 p.K0];

        [time,U,D,V,B,T] = fullmodsim(p,initial_conditions,tf);
    
        time_s = time(2:end);
        AUC_U(initial_cancer,initial_immune) = sum(U(1:end-1).*[time_s-time(1:end-1)]);
        U_vec(initial_cancer,initial_immune) = U(end);
    
    end
    initial_cancer
end



figure
h = heatmap(real(U_vec),'GridVisible','off','CellLabelColor','none')
colm = cbrewer2('PuBu');
colormap(colm);
xlabel('V_0 (log)')
ylabel('U_0 (log)')
set(gca,'FontSize',18)
