function Y = cCDF_Y_tk(K_t,rho_t,Omega_t,g_0,PL_I2I,PL_I2II,Esp,tt,kk,y)

    chi0_t = @(x)   (x)^2/(x^2-(x-1)^2);
    chi1_t = @(x) (x-1)^2/(x^2-(x-1)^2);
    
    A = chi0_t(kk)*meijerGtol(1-2/Esp,[],0,-2/Esp,y/(g_0*PL_I2II(tt))*(  (kk)/K_t(tt))^Esp);
    B = chi1_t(kk)*meijerGtol(1-2/Esp,[],0,-2/Esp,y/(g_0*PL_I2II(tt))*((kk-1)/K_t(tt))^Esp);
    cCDF_Y0_tk = (1-rho_t(tt))*(2/Esp)*(A-B);
    
    cCDF_Y1_tk = 0;
    for tau = 1:(tt-1)
        barXi = prod(Omega_t((tau+1):tt).*PL_I2I((tau+1):tt));
        
        A = chi0_t(kk)*meijerGtol(1-2/Esp,[],[ones(1,tt-tau),0],-2/Esp,...
            (  (kk)/K_t(tt))^Esp/(g_0*PL_I2II(tt))*y/barXi);
        B = chi1_t(kk)*meijerGtol(1-2/Esp,[],[ones(1,tt-tau),0],-2/Esp,...
            ((kk-1)/K_t(tt))^Esp/(g_0*PL_I2II(tt))*y/barXi);
        
%         if (A-B < 0) || isnan(A-B) || (A-B >= 1e3)
%             a = 1;
%         end
        cCDF_Y1_tk = cCDF_Y1_tk + (2/Esp)*(1-rho_t(tau))*prod(rho_t((tau+1):tt))*(A-B);
    end
    Y = cCDF_Y0_tk + cCDF_Y1_tk;
end