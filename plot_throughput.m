% Thanh Luan Nguyen, Tri Nhu Do, and Georges Kaddoum, 
% "Performance Analysis of Multi-user NOMA Wireless-Powered mMTC Networks:
% A Stochastic Geometry Approach," 
% IEEE Transactions on Communications, Oct., 2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blue1  =[0, 0.4470, 0.7410];
orange1=[0.8500, 0.3250, 0.0980];
yellow1=[0.9290, 0.6940, 0.1250];
purple1=[0.4940, 0.1840, 0.5560];
green1 =[0.4660, 0.6740, 0.1880];
cyan1  =[0.3010, 0.7450, 0.9330];
red1   =[0.6350, 0.0780, 0.1840];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ========= XQoM ===========%
ST_TQoM_SIM = (1-OPe2e_BTEH_MTCD_I_SIM)*R_M_TQoMS/(M-1);
ST_TQoM_ANA = (1-OPe2e_BTEH_MTCD_I_ANA)*R_M_TQoMS/(M-1);
ST_PQoM_SIM = (1-OPe2e_BPEH_MTCD_I_SIM)*R_M_PQoMS/(M-1);
ST_PQoM_ANA = (1-OPe2e_BPEH_MTCD_I_ANA)*R_M_PQoMS/(M-1);
for tt = 1:(M-1)
    ST_TQoM_SIM = ST_TQoM_SIM + (1-OPe2e_TQoM_MTCD_II_SIM(tt,:))*R_t_TQoMS/(M-1);
    ST_TQoM_ANA = ST_TQoM_ANA + (1-OPe2e_TQoM_MTCD_II_ANA(tt,:))*R_t_TQoMS/(M-1);
    ST_PQoM_SIM = ST_PQoM_SIM + (1-OPe2e_PQoM_MTCD_II_SIM(tt,:))*R_t_PQoMS/(M-1);
    ST_PQoM_ANA = ST_PQoM_ANA + (1-OPe2e_PQoM_MTCD_II_ANA(tt,:))*R_t_PQoMS/(M-1);
end
% ========= XCoM ===========%
ST_TCoM_SIM = (1-OPe2e_BTEH_MTCD_I_SIM)*R_M_TCoMS/(M-1);
ST_TCoM_ANA = (1-OPe2e_BTEH_MTCD_I_ANA)*R_M_TCoMS/(M-1);
ST_PCoM_SIM = (1-OPe2e_BPEH_MTCD_I_SIM)*R_M_PCoMS/(M-1);
ST_PCoM_ANA = (1-OPe2e_BPEH_MTCD_I_ANA)*R_M_PCoMS/(M-1);
for tt = 1:(M-1)
    for kk = 1:K_t(tt)
        ST_TCoM_SIM = ST_TCoM_SIM + (1-squeeze(OPe2e_TCoM_MTCD_II_SIM(tt,kk,:))')*R_tk_TCoMS(tt,kk)/(M-1);
        ST_TCoM_ANA = ST_TCoM_ANA + (1-squeeze(OPe2e_TCoM_MTCD_II_ANA(tt,kk,:))')*R_tk_TCoMS(tt,kk)/(M-1);
        ST_PCoM_SIM = ST_PCoM_SIM + (1-squeeze(OPe2e_PCoM_MTCD_II_SIM(tt,kk,:))')*R_tk_PCoMS(tt,kk)/(M-1);
        ST_PCoM_ANA = ST_PCoM_ANA + (1-squeeze(OPe2e_PCoM_MTCD_II_ANA(tt,kk,:))')*R_tk_PCoMS(tt,kk)/(M-1);
    end
end
%% ========================================================================
% Save for CoM/QoM w/o EH
% save('ST_TQoM_NEH_ANA.mat','ST_TQoM_ANA');
% save('ST_TCoM_NEH_ANA.mat','ST_TCoM_ANA');
% save('ST_PQoM_NEH_ANA.mat','ST_PQoM_ANA');
% save('ST_PCoM_NEH_ANA.mat','ST_PCoM_ANA');
%
% save('ST_TQoM_NEH_SIM.mat','ST_TQoM_SIM');
% save('ST_TCoM_NEH_SIM.mat','ST_TCoM_SIM');
% save('ST_PQoM_NEH_SIM.mat','ST_PQoM_SIM');
% save('ST_PCoM_NEH_SIM.mat','ST_PCoM_SIM');
%% ========================================================================
% Save for CoM/QoM w/o EH
% save('ST_CNRR_ANA.mat','ST_TQoM_ANA');
% save('ST_CNRR_SIM.mat','ST_TQoM_SIM');
%% ========================================================================
% Save for CoM/QoM w/o EH
% save('ST_TQoM_ANA.mat','ST_TQoM_ANA');
% save('ST_TCoM_ANA.mat','ST_TCoM_ANA');
% save('ST_PQoM_ANA.mat','ST_PQoM_ANA');
% save('ST_PCoM_ANA.mat','ST_PCoM_ANA');
%
% save('ST_TQoM_SIM.mat','ST_TQoM_SIM');
% save('ST_TCoM_SIM.mat','ST_TCoM_SIM');
% save('ST_PQoM_SIM.mat','ST_PQoM_SIM');
% save('ST_PCoM_SIM.mat','ST_PCoM_SIM');
%% ========================================================================
% Load for CoM/QoM w/o EH
ST_TQoM_NEH_ANA = cell2mat(struct2cell(load('ST_TQoM_NEH_ANA.mat')));
ST_TCoM_NEH_ANA = cell2mat(struct2cell(load('ST_TCoM_NEH_ANA.mat')));
ST_PQoM_NEH_ANA = cell2mat(struct2cell(load('ST_PQoM_NEH_ANA.mat')));
ST_PCoM_NEH_ANA = cell2mat(struct2cell(load('ST_PCoM_NEH_ANA.mat')));

ST_TQoM_NEH_SIM = cell2mat(struct2cell(load('ST_TQoM_NEH_SIM.mat')));
ST_TCoM_NEH_SIM = cell2mat(struct2cell(load('ST_TCoM_NEH_SIM.mat')));
ST_PQoM_NEH_SIM = cell2mat(struct2cell(load('ST_PQoM_NEH_SIM.mat')));
ST_PCoM_NEH_SIM = cell2mat(struct2cell(load('ST_PCoM_NEH_SIM.mat')));
% %% ========================================================================
% % % Load for CoM/QoM w/o EH
ST_CNRR_ANA = cell2mat(struct2cell(load('ST_CNRR_ANA.mat')));
ST_CNRR_SIM = cell2mat(struct2cell(load('ST_CNRR_SIM.mat')));
% %% ========================================================================
% % Load for Regenerative Relaying
ST_TQoM_ANA = cell2mat(struct2cell(load('ST_TQoM_ANA.mat')));
ST_TCoM_ANA = cell2mat(struct2cell(load('ST_TCoM_ANA.mat')));
ST_PQoM_ANA = cell2mat(struct2cell(load('ST_PQoM_ANA.mat')));
ST_PCoM_ANA = cell2mat(struct2cell(load('ST_PCoM_ANA.mat')));

ST_TQoM_SIM = cell2mat(struct2cell(load('ST_TQoM_SIM.mat')));
ST_TCoM_SIM = cell2mat(struct2cell(load('ST_TCoM_SIM.mat')));
ST_PQoM_SIM = cell2mat(struct2cell(load('ST_PQoM_SIM.mat')));
ST_PCoM_SIM = cell2mat(struct2cell(load('ST_PCoM_SIM.mat')));
%% ========================================================================
figure;
tqom_ana = plot(P0dBm_ANA,ST_TQoM_ANA,'-','color',   red1,'linewidth',1); hold on;
tcom_ana = plot(P0dBm_ANA,ST_TCoM_ANA,'-','color',  blue1,'linewidth',1); hold on;
pqom_ana = plot(P0dBm_ANA,ST_PQoM_ANA,'-','color', green1,'linewidth',1); hold on;
pcom_ana = plot(P0dBm_ANA,ST_PCoM_ANA,'-','color',purple1,'linewidth',1); hold on;

tqom_neh_ana = plot(P0dBm_ANA,ST_TQoM_NEH_ANA,'--','color',   red1,'linewidth',1); hold on;
plot(P0dBm_ANA,ST_TQoM_NEH_ANA(end)*ones(size(P0dBm_ANA)),'--','color',   red1); hold on;
tcom_neh_ana = plot(P0dBm_ANA,ST_TCoM_NEH_ANA,'--','color',  blue1,'linewidth',1); hold on;
plot(P0dBm_ANA,ST_TCoM_NEH_ANA(end)*ones(size(P0dBm_ANA)),'--','color',   blue1); hold on;

cnrr_ana = plot(P0dBm_ANA,ST_CNRR_ANA,'-k','linewidth',1); hold on;
plot(P0dBm_ANA,ST_CNRR_ANA(end)*ones(size(P0dBm_ANA)),'--k'); hold on;

SIM = 1:2:length(P0dBm_SIM);
tqom_sim = plot(P0dBm_SIM(SIM),ST_TQoM_SIM(SIM),'o','color','k',...
    'markerfacecolor','None','markersize',5,'linewidth',0.5); hold on;
tcom_sim = plot(P0dBm_SIM(SIM),ST_TCoM_SIM(SIM),'o','color','k',...
    'markerfacecolor','None','markersize',5,'linewidth',0.5); hold on;
pqom_sim = plot(P0dBm_SIM(SIM),ST_PQoM_SIM(SIM),'o','color','k',...
    'markerfacecolor','None','markersize',5,'linewidth',0.5); hold on;
pcom_sim = plot(P0dBm_SIM(SIM),ST_PCoM_SIM(SIM),'o','color','k',...
    'markerfacecolor','None','markersize',5,'linewidth',0.5); hold on;

tqom_sim = plot(P0dBm_SIM(SIM),ST_TQoM_NEH_SIM(SIM),'o','color','k',...
    'markerfacecolor','None','markersize',5,'linewidth',0.5); hold on;
tcom_sim = plot(P0dBm_SIM(SIM),ST_TCoM_NEH_SIM(SIM),'o','color','k',...
    'markerfacecolor','None','markersize',5,'linewidth',0.5); hold on;

plot(P0dBm_SIM(SIM),ST_CNRR_SIM(SIM),'o','color','k',...
    'markerfacecolor','None','markersize',5,'linewidth',0.5); hold on;

leg = legend([tqom_ana tcom_ana pqom_ana pcom_ana ...
    tqom_neh_ana tcom_neh_ana cnrr_ana tqom_sim],...
    {'TQoM (ana.)','TCoM (ana.)',...
    'PQoM (ana.)','PCoM (ana.)',...
    'QoM w/o EH (ana.)',...
    'CoM w/o EH (ana.)',...
    'CNRR','sim'},...
    'NumColumns',3,...
    'Location','best');
xlabel('$P_0$ [dBm]','Interpreter','LaTex');
ylabel('Throughput [bits/s/Hz]');
axis([-50 50 0 1.2]);
% set(gcf,'Position',[100 100 340 330]);
% set(gca, 'LooseInset', get(gca, 'TightInset'));
% leg.ItemTokenSize = [15,18];