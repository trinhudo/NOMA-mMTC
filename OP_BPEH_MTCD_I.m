function out = OP_BPEH_MTCD_I(M,PL_I2I,r_t,PA_QoMS,CPA_QoMS,...
    rho_t,eta_t,beta_t,g_0,lambda_t,R_M_QoMS,OPT)
 
    iota_t = exp(-pi*lambda_t.*r_t.^2); % - Null Probability Observed at DI_0 --> DI_(M-1)
    biota_t= 1 - iota_t;
    Omega_t= beta_t.*eta_t;
    brho_t = 1 - rho_t;
    
    % Outage Probability at Type-I MTCD: (DI_1 --> DI_M)
    SP_MTCDI_PQoM = zeros(1,M);
    for t = 1:(M-1)
        a_1 = 2^( (M-1)*R_M_QoMS )-1;
        a_2 = a_1/( 1-beta_t(t+1) );
        a_3 = a_1/( 1-(a_1+1)*(1-PA_QoMS(2)) );
        a_4 = a_3/( 1-beta_t(t+1) );
        
        switch OPT
            case 'ANA'
                CF_X = @(x) CCDF_X_ANA(rho_t,g_0,PL_I2I,Omega_t,t+1,x);
            case 'ASY'
                CF_X = @(x) CCDF_X_ASY(rho_t,g_0,PL_I2I,Omega_t,t+1,x);
            case 'FLO'
                CF_X = @(x) CCDF_X_FLO(rho_t,t+1,x);
            otherwise
        end
        
        SP_MTCDI_PQoM(t+1) = (brho_t(t+1))*iota_t(t)*CF_X(a_1)...
                            + (rho_t(t+1))*iota_t(t)*CF_X(a_2);
        if (a_1 < PA_QoMS(2)/CPA_QoMS(1))
            SP_MTCDI_PQoM(t+1) = SP_MTCDI_PQoM(t+1)...
                + (brho_t(t+1))*biota_t(t)*CF_X(a_3)...
                +  (rho_t(t+1))*biota_t(t)*CF_X(a_4);
        end
    end
    out = 1 - SP_MTCDI_PQoM;
end