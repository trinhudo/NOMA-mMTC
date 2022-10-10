% tt = 1;
% A = varphi_t(tt,:);
% A(isQoMS(tt,:)==0) = [];
%     
% p = fitdist((A'),'Burr');
% mu_t(tt) = p.alpha;
% theta_t(tt) = p.c;
% m_t(tt) = p.k;
% 
% B = random('Burr',mu_t(tt),theta_t(tt),m_t(tt),[1,nbits]);
% ecdf(A); hold on;
% ecdf(B);
% set(gca,'XScale','log');
% axis([1e-5 1e5 0 1]);

% tt = 3;
% A = d_t(tt,N_t >= 1);
% ksdensity(A); hold on;
% 
% fd = @(r) (2*pi*lambda_t(tt)*r)/(1-exp(-pi*lambda_t(tt)*r_t(tt)^2))...
%     .* exp(-pi*lambda_t(tt).*r.^2);
% fplot(fd,[0,20]);

tt = 2;
kk = 1;

varphi = squeeze(varphi_tk(tt,kk,:))';
isExist= squeeze(isCoMS(tt,kk,:))';

sinr18 = PA_CoMS(tt,K_t(tt)+1)*P_BPEH_t(tt,:).*varphi...
    ./(CPA_CoMS(tt,1)*P_BPEH_t(tt,:).*varphi+noiseVar);
rate18 = tauBPEH.*log2( 1+sinr18 );

SPEvent = (rate18 >= R_M_CoMS);

qqs = (kk:K_t(tt))+1;
a_1 = 2^( (M-1)*R_M_CoMS )-1;
b_1 = 2.^( (M-1)*R_tk_CoMS(tt,kk:K_t(tt)) )-1;
b_3 = max( [ a_1/(PA_CoMS(tt, K_t(tt)+1) - a_1*CPA_CoMS(tt,1)),...
            b_1./(PA_CoMS(tt,kk:K_t(tt)) - b_1.*CPA_CoMS(tt,qqs))] );
        
SPEvent21= (P_BPEH_t(tt,:).*varphi/noiseVar >= b_3);
OP1(tt,kk) = 1 - mean(SPEvent.*SPEvent21.*isExist)/mean(isExist)

chi0_t = @(x)   (x)^2/(x^2-(x-1)^2);
chi1_t = @(x) (x-1)^2/(x^2-(x-1)^2);
A = chi0_t(kk)*meijerGtol(1-2/pathlosExp,[],0,-2/pathlosExp,b_3/(g_0*PL_I2II(tt))*(  (kk)/K_t(tt))^pathlosExp);
B = chi1_t(kk)*meijerGtol(1-2/pathlosExp,[],0,-2/pathlosExp,b_3/(g_0*PL_I2II(tt))*((kk-1)/K_t(tt))^pathlosExp);
Y0_1 = (1-rho_t(tt))*(2/pathlosExp)*(A-B);

Omega_t = beta_t.*eta_t;
Y1_1 = 0;
for tau = 1:(tt-1)
    barXi = prod(Omega_t((tau+1):tt).*PL_I2I((tau+1):tt));

    A = chi0_t(kk)*meijerGtol(1-2/pathlosExp,[],[ones(1,tt-tau),0],-2/pathlosExp,(  (kk)/K_t(tt))^pathlosExp/(g_0*PL_I2II(tt))*b_3/barXi);
    B = chi1_t(kk)*meijerGtol(1-2/pathlosExp,[],[ones(1,tt-tau),0],-2/pathlosExp,((kk-1)/K_t(tt))^pathlosExp/(g_0*PL_I2II(tt))*b_3/barXi);

    Y1_1 = Y1_1 + (2/pathlosExp)*(1-rho_t(tau))*prod(rho_t((tau+1):tt))*(A-B);
end
OP2 = 1-(Y0_1+Y1_1)
% 
% C = squeeze(varphi_tk(tt,kk,:))';
% isCoMSind = squeeze(isCoMS(tt,kk,:))';
% C(~logical(isCoMSind)) = [];
% 
% figure;
% ecdf(100*C); hold on;
% plot(yy,cCDF_Y0_tk(1:length(yy)),'--');
% set(gca,'XScale','log');
% set(gca,'YScale','log');