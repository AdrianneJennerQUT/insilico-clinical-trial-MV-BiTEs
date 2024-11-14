function [time, U, D, B, K] = modelsimulator_ODE_case3(p,tf,initialconds);

tspan = [0 tf];

sol =ode45(@odefun, tspan, initialconds);

time = sol.x;
U = sol.y(1,:);
D = sol.y(2,:);
B = sol.y(3,:);
K = sol.y(4,:);

function dydt = odefun(t,y,Z)

  U = y(1);
  D = y(2);
  B = y(3);
  K = y(4);

  T = U+D+K;

  dU = p.r*U-p.k*U*K/T*(p.eps+B/(B+p.gam));
  dD = p.k*U*K/T*(p.eps+B/(B+p.gam));
  dB = -p.d_B*B;
  dK = -p.d_K*K;
  
  dydt = [dU; dD; dB; dK];

end

end