function [time, U_mat, D_mat, B_mat, K_mat, Lysis_mod xmultinonlin] = fit_ODE_model_BiTEs(p,tf,BiTE_conc,Lysis)

param_guess = [p.k,p.eps,p.gam];

lb = [0 0 0 ];
ub = [10 10 1e5];

problem = createOptimProblem('lsqnonlin','x0',param_guess,'objective',@residualsfunction,...
    'lb',lb,'ub',ub');

figure
ms = MultiStart('PlotFcns',@gsplotbestf);
[xmultinonlin,errormultinonlin] = run(ms,problem,100);

p.k = xmultinonlin(1);
p.eps = xmultinonlin(2);
p.gam = xmultinonlin(3);


[time, U, D, B, K] = modelsimulator_ODE_case3(p,tf,initialconds);

    initialconds = [1e5 0 BiTE_conc(1) 2.5*1e5];
    [time, U, D, B, K] = modelsimulator_ODE_case3(p,tf,initialconds);
    Lysis_mod(1) = D(end)/(D(end)+U(end));
    U_mat(1,:) = interp1(time,U,linspace(0,tf,100));
    D_mat(1,:) = interp1(time,D,linspace(0,tf,100));
    B_mat(1,:) = interp1(time,B,linspace(0,tf,100));
    K_mat(1,:) = interp1(time,K,linspace(0,tf,100));

    initialconds = [1e5 0 BiTE_conc(2) 2.5*1e5];
    [time, U, D, B, K] = modelsimulator_ODE_case3(p,tf,initialconds);
    Lysis_mod(2) = D(end)/(D(end)+U(end));
    U_mat(2,:) = interp1(time,U,linspace(0,tf,100));
    D_mat(2,:) = interp1(time,D,linspace(0,tf,100));
    B_mat(2,:) = interp1(time,B,linspace(0,tf,100));
    K_mat(2,:) = interp1(time,K,linspace(0,tf,100));

    initialconds = [1e5 0 BiTE_conc(3) 2.5*1e5];
    [time, U, D, B, K] = modelsimulator_ODE_case3(p,tf,initialconds);
    Lysis_mod(3) = D(end)/(D(end)+U(end));
    U_mat(3,:) = interp1(time,U,linspace(0,tf,100));
    D_mat(3,:) = interp1(time,D,linspace(0,tf,100));
    B_mat(3,:) = interp1(time,B,linspace(0,tf,100));
    K_mat(3,:) = interp1(time,K,linspace(0,tf,100));

    initialconds = [1e5 0 BiTE_conc(4) 2.5*1e5];
    [time, U, D, B, K] = modelsimulator_ODE_case3(p,tf,initialconds);
    Lysis_mod(4) = D(end)/(D(end)+U(end));
    U_mat(4,:) = interp1(time,U,linspace(0,tf,100));
    D_mat(4,:) = interp1(time,D,linspace(0,tf,100));
    B_mat(4,:) = interp1(time,B,linspace(0,tf,100));
    K_mat(4,:) = interp1(time,K,linspace(0,tf,100));

    initialconds = [1e5 0 BiTE_conc(5) 2.5*1e5];
    [time, U, D, B, K] = modelsimulator_ODE_case3(p,tf,initialconds);
    Lysis_mod(5) = D(end)/(D(end)+U(end));
    U_mat(5,:) = interp1(time,U,linspace(0,tf,100));
    D_mat(5,:) = interp1(time,D,linspace(0,tf,100));
    B_mat(5,:) = interp1(time,B,linspace(0,tf,100));
    K_mat(5,:) = interp1(time,K,linspace(0,tf,100));

    initialconds = [1e5 0 BiTE_conc(6) 2.5*1e5];
    [time, U, D, B, K] = modelsimulator_ODE_case3(p,tf,initialconds);
    Lysis_mod(6) = D(end)/(D(end)+U(end));
    U_mat(6,:) = interp1(time,U,linspace(0,tf,100));
    D_mat(6,:) = interp1(time,D,linspace(0,tf,100));
    B_mat(6,:) = interp1(time,B,linspace(0,tf,100));
    K_mat(6,:) = interp1(time,K,linspace(0,tf,100));
%-----------------------------------------------------------------------
function  val = residualsfunction(param)

    p.k = param(1);
    p.eps = param(2);
    p.gam = param(3);
   
    initialconds = [1e5 0 BiTE_conc(1) 2.5*1e5];
    [time, U, D, B, K] = modelsimulator_ODE_case3(p,tf,initialconds);
    Lysis_mod(1) = D(end)/(D(end)+U(end));

    initialconds = [1e5 0 BiTE_conc(2) 2.5*1e5];
    [time, U, D, B, K] = modelsimulator_ODE_case3(p,tf,initialconds);
    Lysis_mod(2) = D(end)/(D(end)+U(end));

    initialconds = [1e5 0 BiTE_conc(3) 2.5*1e5];
    [time, U, D, B, K] = modelsimulator_ODE_case3(p,tf,initialconds);
    Lysis_mod(3) = D(end)/(D(end)+U(end));

    initialconds = [1e5 0 BiTE_conc(4) 2.5*1e5];
    [time, U, D, B, K] = modelsimulator_ODE_case3(p,tf,initialconds);
    Lysis_mod(4) = D(end)/(D(end)+U(end));
        
    initialconds = [1e5 0 BiTE_conc(5) 2.5*1e5];
    [time, U, D, B, K] = modelsimulator_ODE_case3(p,tf,initialconds);
    Lysis_mod(5) = D(end)/(D(end)+U(end));

    initialconds = [1e5 0 BiTE_conc(6) 2.5*1e5];
    [time, U, D, B, K] = modelsimulator_ODE_case3(p,tf,initialconds);
    Lysis_mod(6) = D(end)/(D(end)+U(end));

    val = [Lysis_mod*100-Lysis];
    %val(isnan(val)) = [];
    sum(abs(val).^2);
        
end


end