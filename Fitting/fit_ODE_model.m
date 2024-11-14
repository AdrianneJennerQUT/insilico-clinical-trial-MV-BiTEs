function [time, U, V,xmultinonlin,ci] = fit_ODE_model(p,tf,initialconds,cell_data,viral_data,cell_time, viral_time)

param_guess = [p.beta, p.eta, p.d_V, p.alpha];

lb = [1e-5 0 0 0];
ub = [100 1e7 10 viral_data(1)];

problem = createOptimProblem('lsqnonlin','x0',param_guess,'objective',@residualsfunction,...
    'lb',lb,'ub',ub');

figure
ms = MultiStart('PlotFcns',@gsplotbestf);
[xmultinonlin,errormultinonlin,exitflag,output] = run(ms,problem,200);

options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','MaxIterations',400,'MaxFunctionEvaluations',5000); % you can add 'display','iter' here to get it to print out what it's doing during the fitting
[xmultinonlin_2,resnorm,residual,exitflag,fit_out,lambda,jacobian] = lsqnonlin(@residualsfunction,xmultinonlin,[lb],[ub]);%,options);%fmincon(@residualsfunction, parameter_guess, [0 0 0 0], [Inf Inf Inf Inf], options);  % can use fmincon/ lsqnonlin
xmultinonlin-xmultinonlin_2;
ci = nlparci(xmultinonlin_2,residual,'jacobian',jacobian);

p.beta = xmultinonlin(1);
p.eta = xmultinonlin(2);
p.d_V = xmultinonlin(3);
p.r = 0.004;% xmultinonlin(4);
p.h = 1;%xmultinonlin(4);
p.alpha =xmultinonlin(4);

[time, U, V] = modelsimulator_ODE_case1(p,tf,initialconds);

%-----------------------------------------------------------------------
function  val = residualsfunction(param)

    p.beta = param(1);
    p.eta = param(2);
    p.d_V = param(3);
    p.h = 1;%param(4);
    p.alpha = param(4);
   
    [time, U, V] = modelsimulator_ODE_case1(p,tf,initialconds);
    
    Usol = interp1(time,U,cell_time);
    Vsol = interp1(time,V,viral_time);

    if isempty(find(Usol>1e10))==0
       Usol  =repmat(1e10, 1, length(cell_time))';
       Vsol = repmat(1e10, 1, length(viral_time));
    elseif sum(isnan(Usol))>0
       Usol  =repmat(1e10, 1, length(cell_time))';
       Vsol = repmat(1e10, 1, length(viral_time));
    elseif sum(isnan(Vsol))>0
       Usol  =repmat(1e10, 1, length(cell_time))';
       Vsol = repmat(1e10, 1, length(viral_time));
    end

    Usol = interp1(time,U,cell_time);
    cell_viab = (Usol)/(initialconds(1))*100;

    Vsol = interp1(time,V,viral_time);
    val_cell_viab = [cell_viab'- cell_data]./max(cell_data);

    val_viral_proj = [Vsol - viral_data./p.alpha]/max(viral_data);
    
    val = [val_cell_viab val_viral_proj];

    %val(isnan(val)) = [];
    sum(abs(val).^2);
        
end


end