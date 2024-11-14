function [time, U, D, V, B, T] = modelsimulator_ODE_fullmod(p,tf,initialconds);

tspan = [0 tf];

sol =ode45(@odefun, tspan, initialconds);

time = sol.x;
U = sol.y(1,:);
D = sol.y(2,:);
V = sol.y(3,:);
B = sol.y(4,:);
T = sol.y(5,:);

function dydt = odefun(t,y,Z)

  U = y(1);
  D = y(2);
  V = y(3);
  B = y(4);
  T = y(5);

  %if V<1
  %  V = 0;
  %end

  N = U+D+T;

  dU = p.r*U*log(p.L/U)-p.beta*U^p.h*V/(U^p.h+p.eta^p.h)-p.k*U*T/N*(p.eps+B/(B+p.gam)); %*log(p.L/U)
  dD = p.beta*U^p.h*V/(U^p.h+p.eta^p.h)+p.k*U*T/N*(p.eps+B/(B+p.gam))-p.d_D*D;
  dV = p.beta*U^p.h*V/(U^p.h+p.eta^p.h)-p.d_V*V;
  dB = p.alpha_2*p.beta*U^p.h*V/(U^p.h+p.eta^p.h)-p.d_B*B;
  dT = p.s*D-p.d_T*T;
  
  dydt = [dU; dD; dV; dB; dT];

end

end