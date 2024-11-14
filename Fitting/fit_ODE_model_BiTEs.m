function [time, U, V, B, xmultinonlin] = fit_ODE_model_BiTEs(p,tf,initialconds,BITE_data,BITE_time)

param_guess = [p.alpha_2 1];

lb = [0 0];
ub = [10 10];

problem = createOptimProblem('lsqnonlin','x0',param_guess,'objective',@residualsfunction,...
    'lb',lb,'ub',ub');

figure
ms = MultiStart('PlotFcns',@gsplotbestf);
[xmultinonlin,errormultinonlin] = run(ms,problem,100);

p.alpha_2 = xmultinonlin(1);
p.d_B = xmultinonlin(2);

[time, U, V, B] = modelsimulator_ODE_case2(p,tf,initialconds);

%-----------------------------------------------------------------------
function  val = residualsfunction(param)

    p.alpha_2 = param(1);
    p.d_B = param(2);

    [time, U, V, B] = modelsimulator_ODE_case2(p,tf,initialconds);
    
    Bsol = interp1(time,B,BITE_time);
    
    val = [Bsol-BITE_data];
    %val(isnan(val)) = [];
    sum(abs(val).^2);
        
end


end