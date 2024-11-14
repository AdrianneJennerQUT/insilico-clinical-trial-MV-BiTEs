function [U_mat, D_mat, V_mat, B_mat, T_mat, U_end, D_end, V_end, B_end, T_end] = model_sim_varying_sched(p,tf,initialconds,totaldosages,daysapart)

        U_mat = [];
        D_mat = [];
        V_mat = [];
        B_mat = [];
        T_mat = [];

if totaldosages ==1
    
        [time, U, D, V, B, T] = modelsimulator_ODE_fullmod_V2(p,tf,initialconds);
        U_mat = U;
        D_mat = D;
        V_mat = V;
        B_mat = B;
        T_mat = T;
else
    initialconds_d = initialconds;
    for i = 1:totaldosages
        [time, U, D, V, B, T] = modelsimulator_ODE_fullmod_V2(p,daysapart,initialconds_d);
        initialconds_d = [U(end),D(end),V(end)+1,B(end),T(end)];
        U_mat = [U_mat U(1:end-1)'];
        D_mat = [D_mat D(1:end-1)'];
        V_mat = [V_mat V(1:end-1)'];
        B_mat = [B_mat B(1:end-1)'];
        T_mat = [T_mat T(1:end-1)'];
    end
    if totaldosages*daysapart<tf
        [time, U, D, V, B, T] = modelsimulator_ODE_fullmod_V2(p,tf-totaldosages*daysapart,initialconds_d);
        U_mat = [U_mat U'];
        D_mat = [D_mat D'];
        V_mat = [V_mat V'];
        B_mat = [B_mat B'];
        T_mat = [T_mat T'];
    end
end

U_end = real(U(end));
D_end = real(D(end));
V_end = real(V(end));
B_end = real(B(end));
T_end = real(T(end));



end

