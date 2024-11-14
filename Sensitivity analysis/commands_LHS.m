%% initialise parameters

%load parameters
p = parameters_invivo();

orig_params(1) =p.beta;
orig_params(2) =p.k;
orig_params(3) =p.eps;
orig_params(4) =p.gamma;
orig_params(5) =p.d_D;
orig_params(6) =p.alpha_B;
orig_params(7) =p.d_V;
orig_params(8) =p.d_B;
orig_params(9) =p.d_T;
orig_params(10) =p.s;
orig_params(11) =p.eta;

tf = 25;

initialconds = [p.U0 0 p.V0 p.B0 p.K0];

LB = 0.01;
UB = 2;

betaLB = orig_params(1)*LB;
betaUB = orig_params(1)*UB;
kLB = orig_params(2)*LB;
kUB = orig_params(2)*UB;
epsLB = orig_params(3)*LB;
epsUB = orig_params(3)*UB;
gamLB = orig_params(4)*LB;
gamUB = orig_params(4)*UB;
d_DLB = orig_params(5)*LB;
d_DUB = orig_params(5)*UB;
alpha_2LB = orig_params(6)*LB;
alpha_2UB = orig_params(6)*UB;
d_VLB = orig_params(7)*LB;
d_VUB = orig_params(7)*UB;
d_BLB = orig_params(8)*LB;
d_BUB = orig_params(8)*UB;
d_TLB = orig_params(9)*LB;
d_TUB = orig_params(9)*UB;
sLB = orig_params(10)*LB;
sUB = orig_params(10)*UB;
etaLB = orig_params(11)*LB;
etaUB = orig_params(11)*UB;

%%

n = 5000;
para_num = 11;

X = lhsdesign(n,para_num); %
    
for i = 1:n
   
    i
    
    p.beta = X(i,1)*(betaUB-betaLB)+betaLB;
    p.k = X(i,2)*(kUB-kLB)+kLB;
    p.eps = X(i,3)*(epsUB-epsLB)+epsLB;
    p.gamma = X(i,4)*(gamUB-gamLB)+gamLB;
    p.d_D = X(i,5)*(d_DUB-d_DLB)+d_DLB;
    p.alpha_B = X(i,6)*(alpha_2UB-alpha_2LB)+alpha_2LB;
    p.d_V = X(i,7)*(d_VUB-d_VLB)+d_VLB;
    p.d_B = X(i,8)*(d_BUB-d_BLB)+d_BLB;
    p.d_T = X(i,9)*(d_TUB-d_TLB)+d_TLB;
    p.s = X(i,10)*(sUB-sLB)+sLB;
    p.eta = X(i,11)*(etaUB-etaLB)+etaLB;    
    
    [time,U,D,V,B,T] = fullmodsim_single_injection(p,initialconds,tf);
    
    U_vec(i) = U(end);
    time_s = time(2:end);
    AUC_U(i) = sum(U(1:end-1).*[time_s-time(1:end-1)]);
    
    T_vec(i) = T(end);
    B_vec(i) = B(end);
    V_vec(i) = V(end);   
    D_vec(i) = D(end);    
    
end

%%

U_vec_real = AUC_U;%real(U_vec);

pcc = corrcoef(U_vec_real,X(:,1));%X(:,1));
beta_pcc = pcc(1,2);

pcc = corrcoef(U_vec_real,X(:,2));
k_pcc = pcc(1,2);

pcc = corrcoef(U_vec_real,X(:,3));
eps_pcc = pcc(1,2);

pcc = corrcoef(U_vec_real,X(:,4));
gam_pcc = pcc(1,2);

pcc = corrcoef(U_vec_real,X(:,5));
d_D_pcc = pcc(1,2);

pcc = corrcoef(U_vec_real,X(:,6));
alpha_2_pcc = pcc(1,2);

pcc = corrcoef(U_vec_real,X(:,7));
d_V_pcc = pcc(1,2);

pcc = corrcoef(U_vec_real,X(:,8));
d_B_pcc = pcc(1,2);

pcc = corrcoef(U_vec_real,X(:,9));
d_T_pcc = pcc(1,2);

pcc = corrcoef(U_vec_real,X(:,10));
s_pcc = pcc(1,2);

pcc = corrcoef(U_vec_real,X(:,11));
eta_pcc = pcc(1,2);


%% PCC plot for U

figure
bar([beta_pcc, k_pcc, eps_pcc, gam_pcc, d_D_pcc,...
    alpha_2_pcc, d_V_pcc, d_B_pcc, d_T_pcc, s_pcc, eta_pcc],'FaceColor',[171,221,164]/255,'EdgeColor',[102,194,165]/255)
ylabel('PCC')
set(gca,'FontSize',18)
set(gca,'xtick',[1:11],'xticklabels',...
    {'\beta','k','\epsilon','\gamma','d_D','\alpha_B','d_V','d_B','d_T','s','\eta'})

%%










