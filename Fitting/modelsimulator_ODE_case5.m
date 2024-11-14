function [time, modsol] = modelsimulator_ODE_case5(p,tf,initialconds);

tspan = [0 1];
sol1 =ode45(@odefun, tspan, initialconds);

tspan = [1 2];
sol2 =ode45(@odefun, tspan, sol1.y(:,end)'+[0 0 0 0]);

tspan = [2 3];
sol3 =ode45(@odefun, tspan, sol2.y(:,end)'+[0 0 0 0]);

tspan = [3 4];
sol4 =ode45(@odefun, tspan, sol3.y(:,end)'+[0 0 0 0]);

tspan = [4 tf];
sol5 =ode45(@odefun, tspan, sol4.y(:,end)'+[0 0 0 0]);

time = [sol1.x(1:end-1) sol2.x(1:end-1) sol3.x(1:end-1) sol4.x(1:end-1) sol5.x];

modsol = [sol1.y(:,1:end-1)'; sol2.y(:,1:end-1)'; sol3.y(:,1:end-1)'; sol4.y(:,1:end-1)'; sol5.y'];


function dydt = odefun(t,y,Z)

  U = y(1);
  D = y(2);
  B = y(3);
  T = y(4);
    
  N = U+D+T;

  dU = p.r*log(p.L/U)*U-p.k*U*T/N*(p.eps+B/(B+p.gamma));
  dD = p.k*U*T/N*(p.eps+B/(B+p.gamma))-p.d_D*D;
  dB = -p.d_B*B;
  dT = p.s*D-p.d_T*T;
  
  dydt = [dU;dD;dB;dT];

end

end