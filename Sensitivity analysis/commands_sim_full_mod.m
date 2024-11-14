
%%simulate full model

%load parameters
p = parameters_invivo();

%set initial condition
initial_conditions = [p.U0, 0, p.V0, p.B0, p.K0];

%time
tf = 30;

% solve model
[time,U,D,V,B,T] = fullmodsim(p,initial_conditions,tf);

%plot
figure
subplot(2,2,1)
hold on 
plot(time,(U+D)/1e6)
ylim([0 1200])
subplot(2,2,2)
hold on 
plot(time,V)
subplot(2,2,3)
hold on 
plot(time,B)
subplot(2,2,4)
hold on 
plot(time,T)

dose_V = 1;
dose_B = 0;
[time,U,D,V,B,T] = fullmodsim_dosagevary(p,initial_conditions,tf,dose_V,dose_B);

subplot(2,2,1)
hold on 
plot(time,(U+D)/1e6)
