function [time,U,D,V,B,T] = fullmodsim_single_injection(p,initial_conditions,tf);


[sol_1] = ode45(@(t,z)mod(t,z,p),[0 tf],initial_conditions);

time = [sol_1.x];
U = [sol_1.y(1,:)];
D = [sol_1.y(2,:)];
V = [sol_1.y(3,:)];
B = [sol_1.y(4,:)];
T = [sol_1.y(5,:)];

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