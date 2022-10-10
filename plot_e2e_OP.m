% Thanh Luan Nguyen, Tri Nhu Do, and Georges Kaddoum, 
% "Performance Analysis of Multi-user NOMA Wireless-Powered mMTC Networks:
% A Stochastic Geometry Approach," 
% IEEE Transactions on Communications, Oct., 2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blue1  =[0, 0.4470, 0.7410];
orange1=[0.8500, 0.3250, 0.0980];
yellow1=[0.9290, 0.6940, 0.1250];
purple1=[0.4940, 0.1840, 0.5560];
green1 =[0.4660, 0.6740, 0.1880];
cyan1  =[0.3010, 0.7450, 0.9330];
red1   =[0.6350, 0.0780, 0.1840];
%% ========================================================================
P0dBm_SIM = -60:5:150;
P0dBm_ANA = -60:5:150;
M = 4;
%% ========================================================================
% Execute for Saving Non-Regenerative Relaying (rho_t = 0)
% if (sum(rho_t==0)==M)
%     save('OPe2e_NEH_MTCD_I_ANA.mat','OPe2e_BTEH_MTCD_I_ANA');
%     save('OPe2e_NEH_MTCD_I_SIM.mat','OPe2e_BTEH_MTCD_I_SIM');
%     save('OPe2e_NEH_MTCD_I_ASY.mat','OPe2e_BTEH_MTCD_I_ASY');
% end
%% ========================================================================
% Execute for Saving Non-Regenerative Relaying (rho_t = 0)
% if (sum(rho_t==0)==M)
%     save('OPe2e_CNRR_MTCD_I_ANA.mat','OPe2e_BTEH_MTCD_I_ANA');
%     save('OPe2e_CNRR_MTCD_I_SIM.mat','OPe2e_BTEH_MTCD_I_SIM');
%     save('OPe2e_CNRR_MTCD_I_ASY.mat','OPe2e_BTEH_MTCD_I_ASY');
% end
%% ========================================================================
% Save for Regenerative Relaying
% save('OPe2e_BTEH_MTCD_I_ANA.mat','OPe2e_BTEH_MTCD_I_ANA');
% save('OPe2e_BTEH_MTCD_I_ASY.mat','OPe2e_BTEH_MTCD_I_ASY');
% save('OPe2e_BTEH_MTCD_I_SIM.mat','OPe2e_BTEH_MTCD_I_SIM');
%
% save('OPe2e_BPEH_MTCD_I_ANA.mat','OPe2e_BPEH_MTCD_I_ANA');
% save('OPe2e_BPEH_MTCD_I_ASY.mat','OPe2e_BPEH_MTCD_I_ASY');
% save('OPe2e_BPEH_MTCD_I_SIM.mat','OPe2e_BPEH_MTCD_I_SIM');
%% ========================================================================
% % Execute for Loading Non-EH with CoM/QoM
OPe2e_NEH_MTCD_I_ANA = cell2mat(struct2cell(load('OPe2e_NEH_MTCD_I_ANA.mat')));
OPe2e_NEH_MTCD_I_SIM = cell2mat(struct2cell(load('OPe2e_NEH_MTCD_I_SIM.mat')));
OPe2e_NEH_MTCD_I_ASY = cell2mat(struct2cell(load('OPe2e_NEH_MTCD_I_ASY.mat')));
%% ========================================================================
% Execute for Loading Non-Regenerative Relaying
OPe2e_CNRR_MTCD_I_ANA = cell2mat(struct2cell(load('OPe2e_CNRR_MTCD_I_ANA.mat')));
OPe2e_CNRR_MTCD_I_SIM = cell2mat(struct2cell(load('OPe2e_CNRR_MTCD_I_SIM.mat')));
OPe2e_CNRR_MTCD_I_ASY = cell2mat(struct2cell(load('OPe2e_CNRR_MTCD_I_ASY.mat')));
%% ========================================================================
% Load for Regenerative Relaying
% OPe2e_BTEH_MTCD_I_ANA = cell2mat(struct2cell(load('OPe2e_BTEH_MTCD_I_ANA.mat')));
% OPe2e_BTEH_MTCD_I_ASY = cell2mat(struct2cell(load('OPe2e_BTEH_MTCD_I_ASY.mat')));
% OPe2e_BTEH_MTCD_I_SIM = cell2mat(struct2cell(load('OPe2e_BTEH_MTCD_I_SIM.mat')));
%
% OPe2e_BPEH_MTCD_I_ANA = cell2mat(struct2cell(load('OPe2e_BPEH_MTCD_I_ANA.mat')));
% OPe2e_BPEH_MTCD_I_ASY = cell2mat(struct2cell(load('OPe2e_BPEH_MTCD_I_ASY.mat')));
% OPe2e_BPEH_MTCD_I_SIM = cell2mat(struct2cell(load('OPe2e_BPEH_MTCD_I_SIM.mat')));
%% ========================================================================
% Execute for Loading Traditional multihop
OPe2e_BTEH_MTCD_I_TRA = cell2mat(struct2cell(load('OPe2e_BTEH_MTCD_I_TRA.mat')));
OPe2e_BPEH_MTCD_I_TRA = cell2mat(struct2cell(load('OPe2e_BPEH_MTCD_I_TRA.mat')));
%% ========================================================================
OPe2e_BTEH_MTCD_I_NMA = cell2mat(struct2cell(load('OPe2e_BTEH_MTCD_I_NMA.mat')));
OPe2e_BPEH_MTCD_I_NMA = cell2mat(struct2cell(load('OPe2e_BPEH_MTCD_I_NMA.mat')));
figure;
%
ANA = (P0dBm_ANA <= 110);
%
bteh_ana = semilogy(P0dBm_ANA(ANA),OPe2e_BTEH_MTCD_I_ANA(ANA),'-','color',red1,'linewidth',1); hold on;
bpeh_ana = semilogy(P0dBm_ANA(ANA),OPe2e_BPEH_MTCD_I_ANA(ANA),'-','color',blue1,'linewidth',1); hold on;
%
neh_ana = semilogy(P0dBm_ANA,OPe2e_NEH_MTCD_I_ANA,'-','color',yellow1,'linewidth',1); hold on;
cnrr_ana = semilogy(P0dBm_ANA,OPe2e_CNRR_MTCD_I_ANA,'-','color',green1,'linewidth',1); hold on;
%
bteh_asy= semilogy(P0dBm_ANA,OPe2e_BTEH_MTCD_I_ASY,'-.','color','k','linewidth',0.5);hold on;
bpeh_asy= semilogy(P0dBm_ANA,OPe2e_BPEH_MTCD_I_ASY,'-.','color','k','linewidth',0.5);hold on;
neh_asy = semilogy(P0dBm_ANA,OPe2e_NEH_MTCD_I_ASY,'-.','color','k','linewidth',0.5); hold on;
cnrr_asy = semilogy(P0dBm_ANA,OPe2e_CNRR_MTCD_I_ASY,'-.','color','k','linewidth',0.5); hold on;
%
bteh_tra = semilogy(P0dBm_ANA,OPe2e_BTEH_MTCD_I_TRA,...
    '-o','color',orange1,'markersize',4); hold on;
bpeh_tra = semilogy(P0dBm_ANA,OPe2e_BPEH_MTCD_I_TRA,...
    ':o','color',orange1,'markersize',4); hold on;
%
bteh_nma = semilogy(P0dBm_ANA(ANA),OPe2e_BTEH_MTCD_I_NMA(ANA),...
    '-o','color',cyan1,'markersize',4); hold on;
bpeh_nma = semilogy(P0dBm_ANA(ANA),OPe2e_BPEH_MTCD_I_NMA(ANA),...
    ':o','color',cyan1,'markersize',4); hold on;
%
SIM = 1:2:length(P0dBm_SIM);
bteh_sim = semilogy(P0dBm_SIM(SIM),OPe2e_BTEH_MTCD_I_SIM(SIM),'o','color','k',...
    'markerfacecolor','None','markersize',4,'linewidth',0.5); hold on;
bpeh_sim =semilogy(P0dBm_SIM(SIM),OPe2e_BPEH_MTCD_I_SIM(SIM),'o','color','k',...
    'markerfacecolor','None','markersize',4,'linewidth',0.5); hold on;
neh_sim =semilogy(P0dBm_SIM(SIM),OPe2e_NEH_MTCD_I_SIM(SIM),'o','color','k',...
    'markerfacecolor','None','markersize',4,'linewidth',0.5); hold on;
cnrr_sim =semilogy(P0dBm_SIM(SIM),OPe2e_CNRR_MTCD_I_SIM(SIM),'o','color','k',...
    'markerfacecolor','None','markersize',4,'linewidth',0.5); hold on;
%
leg = legend([bteh_ana, bpeh_ana, neh_ana, cnrr_ana, bteh_asy,...
    bteh_tra, bteh_nma, bpeh_nma, bteh_sim],...
    {'TQoM/TCoM (ana.)',...
    'PQoM/PCoM (ana.)',...
    'QoM/CoM w/o EH (ana.)',...
    'CNRR (ana)','asymp.',...
    'QoM/CoM w/ TS-based EH',...
    'QoM/CoM w/ PS-based EH',...
    'TS-based EH w/o QoM/CoM',...
    'PS-based EH w/o QoM/CoM','sim.'},...
    'numcolumns',2, 'Location','best');
xlabel('$P_0$ [dBm]','Interpreter','LaTex','FontSize',14);
% xticks([-50 0 50 100 150]);
% xticklabels({'-50' '0' '50' '---' '\infty'});
ylabel('e2e OP','FontSize',14);
axis([-Inf 100 10^(-6) 10^(-0)]);

% set(gca,'Fontsize',10);
% set(gcf,'Position',[100 100 400 330]);
% set(gca, 'LooseInset', get(gca, 'TightInset'));
% leg.ItemTokenSize = [16,18];
% grid on;