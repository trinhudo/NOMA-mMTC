% -------------------------------------------------------------------------
% Cumulative Distribution Function of the Normalized Received Power
%   Observed at the Type-I MTCD DI_1 -> DI_M
% -------------------------------------------------------------------------
function cCDF_Z = cCDF_Z_t(rho_t,g_0,Omega_t,PL_I2I,m_t,theta_t,mu_t,tt,z) % m = 0:M-1
    Z0 = (1-rho_t(tt))/gamma(m_t(tt))...
        * meijerGtol(1-m_t(tt),[],0,[],(z/mu_t(tt)/g_0)^theta_t(tt));
    
    Z1 = 0; 
    for tau = 1:(tt-1)
        barXi = prod(Omega_t((tau+1):tt).*PL_I2I((tau+1):tt));
        
        P = [1-m_t(tt) 1];
        Q = [0 1; ones(tt-tau,1) theta_t(tt)*ones(tt-tau,1)];
        Z1 = Z1 + (1-rho_t(tau)) * prod(rho_t((tau+1):tt))...
            * WSCPfoxH(tt-tau+1,1,1,tt-tau+1,P,Q,(z/barXi/mu_t(tt)/g_0)^theta_t(tt));
    end
    Z1 = Z1/gamma(m_t(tt));
    
    cCDF_Z = Z0 + Z1;
end