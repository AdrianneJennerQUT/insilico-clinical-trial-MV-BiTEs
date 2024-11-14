function [time, U, V] = modelsimulator_ODE_case1_v3(p,tf,initialconds);

tspan = [0:0.1:tf];

[t, y] =ode45(@odefun, tspan, initialconds);

time = t;%sol.x;
U = y(:,1);%sol.y(1,:);
V = y(:,2);%sol.y(2,:);

function dydt = odefun(t,y,Z)

  U = y(1);
  V = y(2);

 % dU = p.r*U*log(p.L/(U))-p.beta*U^p.h*V/(U^p.h+p.eta^p.h);
  dU = p.r*U-p.beta*U*V;
  dV = p.alpha*p.beta*U*V-p.d_V*V;
  
  dydt = [dU; dV];

end

end