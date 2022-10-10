function out = OP_PCoM_MTCD_II(M,K_t,rho_t,eta_t,beta_t,g_0,PL_I2I,PL_I2II,Esp,PA_CoMS,CPA_CoMS,R_M_CoMS,R_tk_CoMS)

    Omega_t = beta_t.*eta_t;
    % Outage Probability at Type-I MTCD: (DI_1 --> DI_M)
    SP_PCoM_MTCD_II = zeros(M-1,max(K_t));
    for tt = 1:(M-1)
        
        a_1 = 2^( (M-1)*R_M_CoMS )-1;
        for kk = 1:K_t(tt)
            qqs = (kk:K_t(tt))+1;
            
            b_1 = 2.^( (M-1)*R_tk_CoMS(tt,kk:K_t(tt)) )-1;
            b_3 = max( [ a_1/(PA_CoMS(tt, K_t(tt)+1) - a_1*CPA_CoMS(tt,1)),...
                        b_1./(PA_CoMS(tt,kk:K_t(tt)) - b_1.*CPA_CoMS(tt,qqs))] );
                    
            if (a_1 < PA_CoMS(tt,K_t(tt)+1)/CPA_CoMS(tt,1))...
            	&& ( sum(b_1 < PA_CoMS(tt,kk:K_t(tt))./CPA_CoMS(tt,qqs)) == K_t(tt)-kk+1 )
                SP_PCoM_MTCD_II(tt,kk) = cCDF_Y_tk(K_t,rho_t,Omega_t,g_0,PL_I2I,PL_I2II,Esp,tt,kk,b_3);
            end
        end
    end
    out = 1 - SP_PCoM_MTCD_II;
end