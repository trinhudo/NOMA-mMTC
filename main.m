% Thanh Luan Nguyen, Tri Nhu Do, and Georges Kaddoum, 
% "Performance Analysis of Multi-user NOMA Wireless-Powered mMTC Networks:
% A Stochastic Geometry Approach," 
% IEEE Transactions on Communications, Oct., 2022

clear all
close all

% ====================== Initialize Network Setting ====================== %
disp('Initializing Network...');
tic;
init_mMTC_NOMA;
toc;
% ======================================================================== %

% ========================== Simulation Results ========================== %
disp('Running Simulation Result...');
tic;
P0dBm_SIM = -60:5:150;
% linspace(-60,150,15);
for iP0 = 1:length(P0dBm_SIM)
    P_0 = 10^(P0dBm_SIM(iP0)/10);

    P_BTEH_t = [P_0*ones(1,nbits); zeros(M-1,nbits)]; % Transmit Power of DI_(1) -> DI_(M-1) for BTEH
    P_BPEH_t = [P_0*ones(1,nbits); zeros(M-1,nbits)]; % Transmit Power of DI_(1) -> DI_(M-1) for BPEH
    for tt = 2:(M-1)
        P_BTEH_t(tt,:) = (M-1)*alpha_t(tt)*eta_t(tt)/(1-alpha_t(tt))*P_BTEH_t(tt-1,:).*phi_t(tt,:).*(isEH(tt,:))...
            + P_0*(~isEH(tt,:));
        P_BPEH_t(tt,:) = beta_t(tt)*eta_t(tt)*P_BPEH_t(tt-1,:).*phi_t(tt,:).*(isEH(tt,:))...
            + P_0*(~isEH(tt,:));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TQoM - MTCD_I
    OP_BTEH_MTCD_I_SIM(:,iP0) = OPS_BTEH_MTCD_I(M,PA_QoMS,CPA_QoMS,P_BTEH_t,phi_t,noiseVar,R_M_TQoMS,isExist,tauBTEH);
    % TQoM - MTCD_II
    OP_TQoM_MTCD_II_SIM(:,iP0) = OPS_TQoM_MTCD_II(M,PA_QoMS,CPA_QoMS,P_BTEH_t,varphi_t,noiseVar,R_M_TQoMS,R_t_TQoMS,tauBTEH,isQoMS);
    % TCoM - MTCD_II
    OP_TCoM_MTCD_II_SIM(:,:,iP0) = OPS_TCoM_MTCD_II(M,K_t,PA_CoMS,CPA_CoMS,P_BTEH_t,varphi_tk,noiseVar,R_M_TCoMS,R_tk_TCoMS,tauBTEH,isCoMS);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PQoM - MTCD_I
    OP_BPEH_MTCD_I_SIM(:,iP0) = OPS_BPEH_MTCD_I(M,PA_QoMS,CPA_QoMS,P_BPEH_t,beta_t,phi_t,noiseVar,R_M_PQoMS,isExist,isEH);
    % PQoM - MTCD_II
    OP_PQoM_MTCD_II_SIM(:,iP0) = OPS_PQoM_MTCD_II(M,PA_QoMS,CPA_QoMS,P_BPEH_t,varphi_t,noiseVar,R_M_PQoMS,R_t_PQoMS,tauBPEH,isQoMS);
    % PCoM - MTCD_II
    OP_PCoM_MTCD_II_SIM(:,:,iP0) = OPS_PCoM_MTCD_II(M,K_t,PA_CoMS,CPA_CoMS,P_BPEH_t,varphi_tk,noiseVar,R_M_PCoMS,R_tk_PCoMS,tauBPEH,isCoMS);
end
toc;

% ========================== Analytical Results ========================== %
disp('Running Analytical Result...');
tic;
P0dBm_ANA = -60:5:150;
% 0:2.5:70;
% 0:5:200;
% linspace(0,250,25);
for iP0 = 1:length(P0dBm_ANA)
    P_0 = 10^(P0dBm_ANA(iP0)/10);
    g_0 = P_0/noiseVar;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TQoM - MTCD_I
    OP_BTEH_MTCD_I_ANA(:,iP0) =...
        OP_BTEH_MTCD_I(M,PL_I2I,r_t,PA_QoMS,CPA_QoMS,rho_t,eta_t,alpha_t,g_0,lambda_t,R_M_TQoMS,'ANA');
    OP_BTEH_MTCD_I_ASY(:,iP0) =...
        OP_BTEH_MTCD_I(M,PL_I2I,r_t,PA_QoMS,CPA_QoMS,rho_t,eta_t,alpha_t,g_0,lambda_t,R_M_TQoMS,'ASY');
    % TQoM - MTCD_II
    OP_TQoM_MTCD_II_ANA(:,iP0) =...
        OP_TQoM_MTCD_II(M,PL_I2I,m_t,theta_t,mu_t,PA_QoMS,CPA_QoMS,rho_t,eta_t,alpha_t,g_0,R_M_TQoMS,R_t_TQoMS);
    % TCoM - MTCD_II
    OP_TCoM_MTCD_II_ANA(:,:,iP0) =...
        OP_TCoM_MTCD_II(M,K_t,rho_t,eta_t,alpha_t,g_0,PL_I2I,PL_I2II,pathlosExp,PA_CoMS,CPA_CoMS,R_M_TCoMS,R_tk_TCoMS);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PQoM - MTCD_I
    OP_BPEH_MTCD_I_ANA(:,iP0) =...
        OP_BPEH_MTCD_I(M,PL_I2I,r_t,PA_QoMS,CPA_QoMS,rho_t,eta_t,beta_t,g_0,lambda_t,R_M_PQoMS,'ANA');
    OP_BPEH_MTCD_I_ASY(:,iP0) =...
        OP_BPEH_MTCD_I(M,PL_I2I,r_t,PA_QoMS,CPA_QoMS,rho_t,eta_t,beta_t,g_0,lambda_t,R_M_PQoMS,'ASY');
    % PQoM - MTCD_II
    OP_PQoM_MTCD_II_ANA(:,iP0) =...
        OP_PQoM_MTCD_II(M,PL_I2I,m_t,theta_t,mu_t,PA_QoMS,CPA_QoMS,rho_t,eta_t,beta_t,g_0,R_M_PQoMS,R_t_PQoMS);
    % PCoM - MTCD_II
    OP_PCoM_MTCD_II_ANA(:,:,iP0) =...
        OP_PCoM_MTCD_II(M,K_t,rho_t,eta_t,beta_t,g_0,PL_I2I,PL_I2II,pathlosExp,PA_CoMS,CPA_CoMS,R_M_PCoMS,R_tk_PCoMS);
end
% ======================================================================== %
toc;
% ========================== The e2e OP ========================== %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OPe2e_BTEH_MTCD_I_SIM = OPe2e_BXEH_MTCD_I(M,OP_BTEH_MTCD_I_SIM);
OPe2e_BPEH_MTCD_I_SIM = OPe2e_BXEH_MTCD_I(M,OP_BPEH_MTCD_I_SIM);

OPe2e_TQoM_MTCD_II_SIM= OPe2e_XQoM_MTCD_II(OP_TQoM_MTCD_II_SIM,M,OP_BTEH_MTCD_I_SIM);
OPe2e_PQoM_MTCD_II_SIM= OPe2e_XQoM_MTCD_II(OP_PQoM_MTCD_II_SIM,M,OP_BPEH_MTCD_I_SIM);

OPe2e_TCoM_MTCD_II_SIM= OPe2e_XCoM_MTCD_II(K_t,OP_TCoM_MTCD_II_SIM,M,OP_BTEH_MTCD_I_SIM);
OPe2e_PCoM_MTCD_II_SIM= OPe2e_XCoM_MTCD_II(K_t,OP_PCoM_MTCD_II_SIM,M,OP_BPEH_MTCD_I_SIM);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OPe2e_BTEH_MTCD_I_ANA = OPe2e_BXEH_MTCD_I(M,OP_BTEH_MTCD_I_ANA);
OPe2e_BTEH_MTCD_I_ASY = OPe2e_BXEH_MTCD_I(M,OP_BTEH_MTCD_I_ASY);

OPe2e_BPEH_MTCD_I_ANA = OPe2e_BXEH_MTCD_I(M,OP_BPEH_MTCD_I_ANA);
OPe2e_BPEH_MTCD_I_ASY = OPe2e_BXEH_MTCD_I(M,OP_BPEH_MTCD_I_ASY);

OPe2e_TQoM_MTCD_II_ANA= OPe2e_XQoM_MTCD_II(OP_TQoM_MTCD_II_ANA,M,OP_BTEH_MTCD_I_ANA);
OPe2e_PQoM_MTCD_II_ANA= OPe2e_XQoM_MTCD_II(OP_PQoM_MTCD_II_ANA,M,OP_BPEH_MTCD_I_ANA);

OPe2e_TCoM_MTCD_II_ANA= OPe2e_XCoM_MTCD_II(K_t,OP_TCoM_MTCD_II_ANA,M,OP_BTEH_MTCD_I_ANA);
OPe2e_PCoM_MTCD_II_ANA= OPe2e_XCoM_MTCD_II(K_t,OP_PCoM_MTCD_II_ANA,M,OP_BPEH_MTCD_I_ANA);
% ======================================================================== %
plot_throughput;
% ======================================================================== %
plot_e2e_OP;