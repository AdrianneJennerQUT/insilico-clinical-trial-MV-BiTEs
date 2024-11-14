function [time, U, D, V, B, T, xmultinonlin] = fit_ODE_model_BiTEs_invivo(p,tf,initialconds,Days_data,Tvol_data)

param_guess = [p.s p.k];

lb = [0 0 ];
ub = [1e2 1e2];

problem = createOptimProblem('lsqnonlin','x0',param_guess,'objective',@residualsfunction,...
    'lb',lb,'ub',ub');

figure
ms = MultiStart('PlotFcns',@gsplotbestf);
[xmultinonlin,errormultinonlin] = run(ms,problem,100);

p.s = xmultinonlin(1);
p.k = xmultinonlin(2);

   [time, U, D, V, B, T] = modelsimulator_ODE_fullmod(p,tf,initialconds);

%-----------------------------------------------------------------------
function  val = residualsfunction(param)

    p.s = param(1);
    p.k = param(2);
   
   [time, U, D, V, B, T] = modelsimulator_ODE_fullmod(p,tf,initialconds);
    
    Usol = interp1(time,U,Days_data)/1e6;
    
    val = [Usol-Tvol_data];
    %val(isnan(val)) = [];
    sum(abs(val).^2);
  
        
end


end