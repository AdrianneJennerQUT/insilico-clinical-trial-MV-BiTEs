function [time, U xmultinonlin] = fit_ODE_model_invivo(p,tf,Days,Tvol_mean)

param_guess = [p.L,p.r];

lb = [0 0];
ub = [1e11 1];

problem = createOptimProblem('lsqnonlin','x0',param_guess,'objective',@residualsfunction,...
    'lb',lb,'ub',ub');

figure
ms = MultiStart('PlotFcns',@gsplotbestf);
[xmultinonlin,errormultinonlin] = run(ms,problem,100);

p.L = xmultinonlin(1);
p.r = xmultinonlin(2);

initialconds = [1e5];
[time, U] = modelsimulator_ODE_case4(p,tf,initialconds);

%-----------------------------------------------------------------------
function  val = residualsfunction(param)

    p.L = param(1);
    p.r = param(2);
   
    initialconds = [1e5];
    [time, U] = modelsimulator_ODE_case4(p,tf,initialconds);
    
    U_mod = interp1(time,U,Days)/1e6;

    val = [U_mod-Tvol_mean];
    %val(isnan(val)) = [];
    sum(abs(val).^2);
        
end


end