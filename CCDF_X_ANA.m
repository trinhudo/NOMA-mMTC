% -------------------------------------------------------------------------
% Cumulative Distribution Function of the Normalized Received Power
%   Observed at the Type-I MTCD DI_2 -> DI_M
% -------------------------------------------------------------------------
function out = CCDF_X_ANA(rho_t,g_0,PL_I2I,Omega_t,t,x)
    %
    t = t-1;
    brho_t = 1-rho_t;
    %
    CCDF_X0_t = brho_t(t) * exp(-x/(g_0*PL_I2I(t+1)));
    %
    CCDF_X1_t = 0;
    for tau = (1):(t-1)
        bar_xi = prod(Omega_t((tau+1):(t)).*PL_I2I((tau+1):(t)));
        inp_val = x/( g_0*PL_I2I(t+1)*bar_xi );
        
%         foo = meijerGtol([],[],[ones(1,t-tau),0],[],inp_val);
        %
        if (inp_val <= 1)
            Fs = @(s) gammaz(1-s).^(t-tau).*gammaz(-s).*(inp_val).^(s); L = -0.5;
        else
            Fs = @(s) gammaz(1+s).^(t-tau).*gammaz(s).*(inp_val).^(-s); L = 0.5;
        end
        foo= integral(Fs, L-1i*100, L+1i*100, 'AbsTol',10^(-12),'RelTol',10^(-12))/(2*pi*1i);
        foo= abs(foo);
        %

        CCDF_X1_t= CCDF_X1_t + brho_t(tau)*prod(rho_t((tau+1):(t))) * foo;
    end
    out = CCDF_X0_t+CCDF_X1_t;
end