function out = OP_PQoM_MTCD_II(M,PL_I2I,m_t,theta_t,mu_t,PA_QoMS,CPA_QoMS,...
        rho_t,eta_t,beta_t,g_0,R_M_QoMS,R_t_QoMS)

    Omega_t = beta_t.*eta_t;
    % Outage Probability at Type-I MTCD: (DI_1 --> DI_M)
    SP_MTCD_II_PQoM = zeros(1,M-1);
    for tt = 1:(M-1)        
        a_1 = 2^( (M-1)*R_M_QoMS )-1;
        b_1 = 2^( (M-1)*R_t_QoMS )-1;
        
        if (a_1 < PA_QoMS(2)/CPA_QoMS(1))
            c_1 = max( a_1/( PA_QoMS(2)-a_1*CPA_QoMS(1) ), b_1/PA_QoMS(1) );
            
            SP_MTCD_II_PQoM(tt) = SP_MTCD_II_PQoM(tt)...
                + cCDF_Z_t(rho_t,g_0,Omega_t,PL_I2I,m_t,theta_t,mu_t,tt,c_1);
        end
    end
    out = 1 - SP_MTCD_II_PQoM;
end