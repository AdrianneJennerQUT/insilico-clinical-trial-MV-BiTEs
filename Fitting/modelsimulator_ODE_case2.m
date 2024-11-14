function [time, U, V, B] = modelsimulator_ODE_case2(p,tf,initialconds);

tspan = [0 tf];

sol =ode45(@odefun, tspan, initialconds);

time = sol.x;
U = sol.y(1,:);
V = sol.y(2,:);
B = sol.y(3,:);

function dydt = odefun(t,y,Z)

  U = y(1);
  V = y(2);
  B = y(3);

  dU = p.r*U-p.beta*U^p.h*V/(U^p.h+p.eta^p.h);
  dV = p.beta*U^p.h*V/(U^p.h+p.eta^p.h)-p.d_V*V;
  dB = p.alpha_2*p.beta*U^p.h*V/(U^p.h+p.eta^p.h)-p.d_B*B;
  
  dydt = [dU; dV;dB];

end

end