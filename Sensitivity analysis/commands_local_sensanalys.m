%% initialise parameters

%load parameters
p = parameters_invivo();

tf = 25;

initialconds = [p.U0 0 p.V0 p.B0 p.K0];


 %% Sensitivity analysis

orig_params(1) =p.beta;
orig_params(2) =p.k;
orig_params(3) =p.eps;
orig_params(4) =p.gamma;
orig_params(5) =p.d_D;
orig_params(6) =p.alpha_B;
orig_params(7) =p.d_V;
orig_params(8) =p.d_B;
orig_params(9) =p.d_T;
orig_params(10) =p.s;
orig_params(11) =p.eta;

[time,U,D,V,B,T] = fullmodsim_single_injection(p,initialconds,tf);

%%
original_day_5 =interp1(time,U,5);
original_day_end = U(end);

timegrid = linspace(0,tf,50);

change = [0.7 0.8 0.9 1 1.1 1.2 1.3];

change = [0.5 0.8 0.9 1 1.1 1.2 1.5];

param2change = 1;
[tumour_mat_end_beta] = sim_sen_ana(p,change,param2change,orig_params,tf,initialconds,timegrid);

param2change = 2;
[tumour_mat_end_k] = sim_sen_ana(p,change,param2change,orig_params,tf,initialconds,timegrid);

 param2change = 3;
[tumour_mat_end_eps] = sim_sen_ana(p,change,param2change,orig_params,tf,initialconds,timegrid);

 param2change = 4;
[tumour_mat_end_gam] = sim_sen_ana(p,change,param2change,orig_params,tf,initialconds,timegrid);

 param2change = 5;
[tumour_mat_end_d_D] = sim_sen_ana(p,change,param2change,orig_params,tf,initialconds,timegrid);

 param2change = 6;
[tumour_mat_end_alpha_2] = sim_sen_ana(p,change,param2change,orig_params,tf,initialconds,timegrid);

 param2change = 7;
[tumour_mat_end_d_V] = sim_sen_ana(p,change,param2change,orig_params,tf,initialconds,timegrid);

 param2change = 8;
[tumour_mat_end_d_B] = sim_sen_ana(p,change,param2change,orig_params,tf,initialconds,timegrid);

 param2change = 9;
[tumour_mat_end_d_T] = sim_sen_ana(p,change,param2change,orig_params,tf,initialconds,timegrid);

 param2change = 10;
[tumour_mat_end_s] = sim_sen_ana(p,change,param2change,orig_params,tf,initialconds,timegrid);

 param2change = 11;
[tumour_mat_end_eta] = sim_sen_ana(p,change,param2change,orig_params,tf,initialconds,timegrid);

%%

matrix_plot = (real([tumour_mat_end_beta;tumour_mat_end_k;tumour_mat_end_eps;...
    tumour_mat_end_gam;tumour_mat_end_d_D;tumour_mat_end_alpha_2;...
    tumour_mat_end_d_V;tumour_mat_end_d_B;tumour_mat_end_d_T;tumour_mat_end_s;tumour_mat_end_eta])-original_day_end)/original_day_end;

figure
h = heatmap(matrix_plot,'GridVisible','off')

col_map = cbrewer2('RdBu');
colormap(flipud(col_map))

XLabels = 1:length(change)*11;
% Convert each number in the array into a string
CustomXLabels = string(XLabels);
% Replace all but the fifth elements by spaces
CustomXLabels(mod(XLabels,100) ~= 0) = " ";

CustomXLabels(4) = '\beta';
CustomXLabels(11) = 'k';
CustomXLabels(18) = '\epsilon';
CustomXLabels(25) = '\gamma';
CustomXLabels(32) = 'd_D';
CustomXLabels(39) = '\alpha_V';
CustomXLabels(46) = 'd_V';
CustomXLabels(53) = 'd_B';
CustomXLabels(60) = 'd_T';
CustomXLabels(67) = 's';
CustomXLabels(74) = '\eta';

% Set the 'XDisplayLabels' property of the heatmap 
% object 'h' to the custom x-axis tick labels
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = CustomXLabels;
set(gca,'FontSize',17)


%% 

p = parameters_invivo();

[time,U,D,V,B,T] = fullmodsim_single_injection(p,initialconds,tf);



