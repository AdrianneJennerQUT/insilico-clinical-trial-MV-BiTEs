%% determine impact of varying V0 and B0 size for fixed dosage protocol


%load parameters
p = parameters_invivo();
p.B0tot = p.B0*5;
p.V0tot = p.V0*5;

tf = 30;


%% vary number of dosages and timing of dosages

[time,U,D,V,B,T] = fullmodsim(p,initial_conditions,tf);

%%

for doses = 1:5
    
    for days_apart = 1:5
   
        p.B0 = p.B0tot/doses;
        p.V0 = p.V0tot/doses;
        
        initial_conditions = [p.U0 0 p.V0 p.B0 p.K0];
        %simulate model
        [time,U,D,V,B,T] = fullmodsim_multiplicity_timing_vary(p,initial_conditions,tf,doses,days_apart);
        
        %save AUC and U(end)
        U_vec(doses,days_apart) = U(end);        
        
    end
    
end
        

%%

figure
h = heatmap(U_vec,'GridVisible','off','CellLabelColor','none')

colm = cbrewer2('GnBu');
colormap(colm);
ylabel('# injections')
xlabel('Days apart')
set(gca,'FontSize',18)
colorbar






