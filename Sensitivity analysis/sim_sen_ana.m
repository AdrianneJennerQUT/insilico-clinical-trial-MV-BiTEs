function [tumour_mat_end] = sim_sen_ana(p,change,param2change,orig_params,tf,initialconds,timegrid)

params = orig_params;
tumour_mat_end = NaN(length(change),length(change)*length(params));

for i = 1:length(change)
    
    mat_change_end = NaN(1,length(change)*length(params));

    
    for j = 1:length(params)
       
            if j == param2change
                    params = orig_params;
                    params(param2change) = orig_params(param2change)*change(i);

                    p.beta = params(1);
                    p.k = params(2);
                    p.eps = params(3);
                    p.gamma = params(4);
                    p.d_D = params(5);
                    p.alpha_v = params(6);
                    p.d_V = params(7);
                    p.d_B = params(8);
                    p.d_T = params(9);
                    p.s = params(10);
                    p.eta = params(11);

                    [time,U,D,V,B,T] = fullmodsim_single_injection(p,initialconds,tf);

                    if isreal(U(end))==0
                        disp('Complex')
                    end
                    
                    mat_change_end((j-1)*7+i)=U(end);
            else
                for k = 1:length(change)

                    params = orig_params;
                    params(j) = orig_params(j)*change(k);
                    params(param2change) = orig_params(param2change)*change(i);

                    p.beta = params(1);
                    p.k = params(2);
                    p.eps = params(3);
                    p.gamma = params(4);
                    p.d_D = params(5);
                    p.alpha_v = params(6);
                    p.d_V = params(7);
                    p.d_B = params(8);
                    p.d_T = params(9);
                    p.s = params(10);
                    p.eta = params(11);

                    [time,U,D,V,B,T] = fullmodsim_single_injection(p,initialconds,tf);

                    if isreal(U(end))==0
                        disp('Complex')
                    end
                    
                    mat_change_end((j-1)*7+k)=U(end);
                                                            
                end
            end
       

    
    end
    tumour_mat_end(i,:) = mat_change_end;
i
end


end