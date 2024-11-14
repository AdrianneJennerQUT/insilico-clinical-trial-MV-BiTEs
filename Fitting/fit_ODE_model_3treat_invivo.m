function [timeB, modB,timeMV, modMV, timeMVB, modMVB, xmultinonlin] = fit_ODE_model_3treat_invivo(p,tf,Tvol_BiTEs,Tvol_MV,Tvol_CD20,time_BiTE,time_MV,time_MVB)

param_guess = [p.s p.d_D p.beta p.d_B p.eps p.alpha_2 p.k];%[p.s,p.d_D,p.alpha_2,p.beta,p.d_B,p.d_V];

lb = [0 0 0 0 0 0 0];% 0 0 0];
ub = [20 10 10 10 10 1e3 20];%1e6 450 200 30];

problem = createOptimProblem('lsqnonlin','x0',param_guess,'objective',@residualsfunction,...
    'lb',lb,'ub',ub');

figure
ms = MultiStart('PlotFcns',@gsplotbestf);
[xmultinonlin,errormultinonlin] = run(ms,problem,2000);



    p.s = xmultinonlin(1);
    p.d_D = xmultinonlin(2);
    p.beta = xmultinonlin(3);
    p.d_B = xmultinonlin(4);
    p.eps = xmultinonlin(5);
    p.alpha_2 = xmultinonlin(6);
    p.k = xmultinonlin(7);


initialconds = [p.U0 0 p.B0 p.K0];
[timeB, modB] = modelsimulator_ODE_case5(p,tf,initialconds); %BiTEs
initialconds = [p.U0 0 p.V0 p.K0];
[timeMV, modMV] = modelsimulator_ODE_case6(p,tf,initialconds); %MV
initialconds = [p.U0 0 p.V0 p.B0 p.K0];
[timeMVB, modMVB] = modelsimulator_ODE_case7(p,tf,initialconds); %MV-BiTEs

%-----------------------------------------------------------------------
function  val = residualsfunction(param)

    p.s = param(1);
    p.d_D = param(2);
    p.beta = param(3);
    p.d_B = param(4);
    p.eps = param(5);
    p.alpha_2 = param(6);
    p.k = param(7);
  
    initialconds = [p.U0 0 p.B0 p.K0];
    [timeB, modB] = modelsimulator_ODE_case5(p,tf,initialconds); %BiTEs
    initialconds = [p.U0 0 p.V0 p.K0];
    [timeMV, modMV] = modelsimulator_ODE_case6(p,tf,initialconds); %MV
    initialconds = [p.U0 0 p.V0 p.B0 p.K0];
    [timeMVB, modMVB] = modelsimulator_ODE_case7(p,tf,initialconds); %MV-BiTEs
    
    UB_mod = interp1(timeB,modB(:,1)',time_BiTE)/1e6;
    UMV_mod = interp1(timeMV,modMV(:,1)',time_MV)/1e6;
    UMVB_mod = interp1(timeMVB,modMVB(:,1)',time_MVB)/1e6;

    val = [UB_mod'-Tvol_BiTEs UMV_mod'-Tvol_MV UMVB_mod'-Tvol_CD20];
    %val = [UMV_mod'-Tvol_MV];

    if isnan(val)==1
        val = repmat(1e10,1,length(val))
    end
   
    sum(abs(val).^2);
        
end


end