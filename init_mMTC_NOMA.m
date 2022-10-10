nbits = 1e6;
d = 1; % Temporary Shifted Index/Convert to MatLab Index
M = 4; % Number of Type-I MTCD
B = 10; % MHz - BW
pathlosExp = 3.67;
carrierFreq = 3; % GHz - Carrier Frequency
noisePowerDensity = -174; % dBm/Hz
noiseVar = 10^(noisePowerDensity/10)*B*1e6; % mW - Noise Variance

% l_q = [Inf 75 50 50]; % meter - Distance Between Type-I MTCDs
l_q = [Inf 200 200 200]; % meter - Distance Between Type-I MTCDs  
GI_Tx = 0; % dBi - Transmit Antenna Gain of Type-I MTCD
GI_Rx = 0; % dBi - Receive Antenna Gain of Type-I MTCD    
PL_I2I= pathloss3GPP_UMi(GI_Rx,GI_Tx,carrierFreq,l_q); % - Per Hop Pathloss between Type-I MTCDs

% r_t = [50 25 25]; % meter - Coverage Radius Per DI
r_t = [100*ones(1,M-1)]; % meter - Coverage Radius Per DI
GII_Tx = 5; % dBi - Transmit Antenna Gain of Type-II MTCD
GII_Rx = 5; % dBi - Receive Antenna Gain of Type-II MTCD
PL_I2II = pathloss3GPP_UMi(GII_Rx,GII_Tx,carrierFreq, r_t);

lambda_t = (1e-3)*ones(1,M-1); % Type-II MTCDs/m^2 - Per Hop Type-II MTCDs Density
% lambda_t = 0*ones(1,M-1); % Type-II MTCDs/m^2 - Per Hop Type-II MTCDs Density
maxRateQoS = 1.5;
RateQoS = 0.5;

% Connectivity-oriented multi-hop scheme (CoM)
K_t = [3,2,1]; % MTCD - maximumActive_MTCDTypeII_CoMS len(K_m) <= M-1
PA_CoMS = zeros(M-1,max(K_t)+1); % - Power Allocation CoMS
CPA_CoMS= zeros(size(PA_CoMS)); % - Cumulative Power Allocation CoMS
for tt = 1:(M-1)
    [PA_CoMS(tt,1:(K_t(tt)+1)),CPA_CoMS(tt,1:(K_t(tt)+1))] = powCoeff(K_t(tt));
end

rho = 0.01;
rho_t = [0,rho,rho,0]; % - Probability Energy Harvesting at [DI_1,...,DI_M]

eta_t = [0,ones(1,M-2),0]; % - Energy Harvesting Efficiency

alpha_t= [0,0.1*ones(1,M-1)]; % - Energy Harvesting Ratio /BTEH
beta_t = [0,0.8*ones(1,M-1)]; % - Energy Harvesting Ratio /BPEH

% Quality-oriented multi-hop scheme (QoM)
[PA_QoMS,CPA_QoMS] = powCoeff(1);
R_M_TQoMS = log2(1+PA_QoMS(2)/CPA_QoMS(1))*(1-min(alpha_t(2:(M-1)).*(rho_t(2:(M-1))>0)))/(M-1)*RateQoS;
R_t_TQoMS = maxRateQoS*RateQoS;
R_M_PQoMS = log2(1+PA_QoMS(2)/CPA_QoMS(1))/(M-1)*RateQoS;
R_t_PQoMS = maxRateQoS*RateQoS;

R_M_TCoMS = Inf*ones(M,1);
R_M_PCoMS = Inf*ones(M,1);
R_tk_TCoMS= Inf*ones(M-1,max(K_t));
R_tk_PCoMS= Inf*ones(M-1,max(K_t));
for tt = 1:(M-1)
    R_M_TCoMS(tt+1)= log2( 1+PA_CoMS(tt,K_t(tt)+1)/CPA_CoMS(tt,1) )*(1-alpha_t(tt+1)*(rho_t(tt+1)>0))/(M-1);
    R_M_PCoMS(tt+1)= log2( 1+PA_CoMS(tt,K_t(tt)+1)/CPA_CoMS(tt,1) )/(M-1);
    for kk = 1:K_t(tt)
    	R_tk_TCoMS(tt,kk) = log2( 1+PA_CoMS(tt,kk)/CPA_CoMS(tt,kk+1) )*(1-alpha_t(tt+1)*(rho_t(tt+1)>0))/(M-1);
    	R_tk_PCoMS(tt,kk) = log2( 1+PA_CoMS(tt,kk)/CPA_CoMS(tt,kk+1) )/(M-1);
    end
end
R_M_TCoMS = min(R_M_TCoMS,[],'all')*RateQoS;
R_tk_TCoMS(R_tk_TCoMS == Inf) = maxRateQoS;
R_tk_TCoMS= R_tk_TCoMS*RateQoS;

R_M_PCoMS = min(R_M_PCoMS,[],'all')*RateQoS;
R_tk_PCoMS(R_tk_PCoMS == Inf) = maxRateQoS;
R_tk_PCoMS= R_tk_PCoMS*RateQoS;

isQoMS = zeros(M-1,nbits);
isCoMS = zeros(M-1,max(K_t),nbits);
d_t = Inf*ones(M-1,nbits); 
d_tk= Inf*ones(M-1,max(K_t),nbits);
for tt = 1:(M-1)
    N_t = poissrnd(lambda_t(tt)*r_t(tt)^2*pi,[1,nbits]); % Type-II MTCD - Number of Type-II MTCD Exist within Proximity of Type-I MTCDs
    isQoMS(tt,:) = (N_t >= 1); % - Identity Variable Specifying Whether There Exists at Least 1 Type-I MTCD
    d_t(tt,N_t >= 1) = arrayfun(@(x) r_t(tt)*min(sqrt(rand(1,N_t(x)))),find(N_t >= 1)); % meter - Distance from Type-I MTCDs to the Paired Type-II MTCDs/QoM
    d_t(tt,N_t == 0) = Inf;
    
    for kk = 1:K_t(tt)
        Rin = r_t(tt)*(kk-1)/K_t(tt);
        Rout= r_t(tt)*(kk)/K_t(tt);
        N_tk = poissrnd(lambda_t(tt)/K_t(tt)*r_t(tt)^2*pi,[1,nbits]) >= 1; % Type-II MTCD - Number of Type-II MTCD Exist within Proximity of Type-I MTCDs
        isCoMS(tt,kk,:) = N_tk;
        N_tk(N_tk == 0) = Inf;
        
        d_tk(tt,kk,:) = sqrt((Rout^2-Rin^2)*rand(1,nbits)+Rin^2).*N_tk; % meter - Distance from Type-I MTCDs to the Paired Type-II MTCDs/CoM
    end
end
isExist = isQoMS;

PL_I2II_QoM = pathloss3GPP_UMi(GII_Rx,GII_Tx,carrierFreq, d_t); % - Per Hop Pathloss between Type-I MTCD and Type-II MTCDs/QoM
PL_I2II_CoM = pathloss3GPP_UMi(GII_Rx,GII_Tx,carrierFreq,d_tk); % - Per Hop Pathloss between Type-I MTCD and Type-II MTCDs/CoM

isEH = rand(M,nbits) < vec2mat(rho_t,nbits); % - Energy Harvesting Decision Maxtrix
tauBTEH = (1-vec2mat(alpha_t.*(rho_t>0),nbits).*isEH)/(M-1); % - the tau coefficient for BTEH
tauBPEH = 1/(M-1); % - the tau coefficient for BPEH
   
phi_t = [zeros(1,nbits); cell2mat( arrayfun(@(x) PL_I2I(x)*exprnd(1,[1,nbits]),(2:M)','UniformOutput',false) )]; % Normalized Received Pow at Type-I MTCD
varphi_tk= PL_I2II_CoM.*exprnd(1,[M-1,max(K_t),nbits]); % Normalized Received Pow
varphi_t = PL_I2II_QoM.*exprnd(1,[M-1,nbits]); % Normalized Received Pow

% AA = zeros(M-1,nbits);
for tt = 1:M-1
    A = varphi_t(tt,:);
    A(isQoMS(tt,:)==0) = [];
    
    if lambda_t(1) > 0
        p = fitdist(A','Burr');
        mu_t(tt) = p.alpha;
        m_t(tt) = p.k;
        theta_t(tt) = p.c;
    else
        mu_t(tt) = 1;
        m_t(tt) = 1;
        theta_t(tt) = 1;
    end
end