function [time, U] = modelsimulator_ODE_case4(p,tf,initialconds);

tspan = [0 tf];

sol =ode45(@odefun, tspan, initialconds);

time = sol.x;
U = sol.y(1,:);

function dydt = odefun(t,y,Z)

  U = y(1);
    
  dU = p.r*log(p.L/U)*U;
  
  dydt = [dU];

end

end