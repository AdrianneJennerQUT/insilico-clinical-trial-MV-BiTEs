function [time,U,D,V,B,T] = fullmodsim(p,initial_conditions,tf);


[sol_1] = ode45(@(t,z)mod(t,z,p),[0 1],initial_conditions);
[sol_2] = ode45(@(t,z)mod(t,z,p),[1 2],sol_1.y(:,end)'+[0 0 p.V0 p.B0 0]);
[sol_3] = ode45(@(t,z)mod(t,z,p),[2 3],sol_2.y(:,end)'+[0 0 p.V0 p.B0 0]);
[sol_4] = ode45(@(t,z)mod(t,z,p),[3 4],sol_3.y(:,end)'+[0 0 p.V0 p.B0 0]);
[sol_5] = ode45(@(t,z)mod(t,z,p),[4 tf],sol_4.y(:,end)'+[0 0 p.V0 p.B0 0]);

time = [sol_1.x sol_2.x sol_3.x sol_4.x sol_5.x];
U = [sol_1.y(1,:) sol_2.y(1,:) sol_3.y(1,:) sol_4.y(1,:) sol_5.y(1,:)];
D = [sol_1.y(2,:) sol_2.y(2,:) sol_3.y(2,:) sol_4.y(2,:) sol_5.y(2,:)];
V = [sol_1.y(3,:) sol_2.y(3,:) sol_3.y(3,:) sol_4.y(3,:) sol_5.y(3,:)];
B = [sol_1.y(4,:) sol_2.y(4,:) sol_3.y(4,:) sol_4.y(4,:) sol_5.y(4,:)];
T = [sol_1.y(5,:) sol_2.y(5,:) sol_3.y(5,:) sol_4.y(5,:) sol_5.y(5,:)];

function dydt = mod(t,y,p)

    U = y(1);
    D = y(2);
    V = y(3);
    B = y(4);
    T = y(5);
    
	
    if V<1e-16
        V = 0;
    end
    
    if B<1e-16
        B = 0;
    end
    
    if T<1e-16
        T = 0;
    end
    
    if D<1e-16
        D = 0;
    end
    
    
    if U<1e-20
       U = 0;
       dU = 0;
       N = U+D+T;
    else
        N = U+D+T;
        dU = p.r*U*log(p.L/U)-p.beta*U*V/(U+p.eta)-p.k*U*T/N*(p.eps+B/(B+p.gamma));
    end
        
    dD = p.beta*U*V/(U+p.eta)+p.k*U*T/N*(p.eps+B/(B+p.gamma))-p.d_D*D;
    dV = p.beta*U*V/(U+p.eta)-p.d_V*V;
    dB = p.alpha_B*p.beta*U*V/(U+p.eta)-p.d_B*B;
    dT = p.s*D-p.d_T*T;
        

    dydt = [dU;dD;dV;dB;dT];

end


end