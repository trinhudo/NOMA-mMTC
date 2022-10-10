% -------------------------------------------------------------------------
% Cumulative Distribution Function of the Normalized Received Power
%   Observed at the Type-I MTCD DI_2 -> DI_M
% -------------------------------------------------------------------------
function out = CCDF_X_ASY(rho_t,g_0,PL_I2I,Omega_t,tt,x)
    % initialize
    tt = tt - 1;
    brho_t = 1 - rho_t;
    syms s
    % main
    gamcoeff = 1/PL_I2I(tt+1);
    A = rho_t(tt);
    sum_tau = 0;
    for tau = 1:(tt-1)
        sum_tau = sum_tau + brho_t(tau)*prod( rho_t((tau+1):tt) );
    end
    A = A - sum_tau;
    %
    B = brho_t(tt);
    sum_tau = 0;
    for tau = 1:(tt-1)
        prod_j = brho_t(tau)*prod( rho_t((tau+1):tt) );
        bar_xi = prod(Omega_t((tau+1):tt).*PL_I2I((tau+1):tt));
        
        sum_n = 0;
        for n = 0:(tt-tau)
            sum_r = 0;
            for r = 0:(tt-tau-n)
                sum_r = sum_r + (-1)^r/factorial(r)...
                    * log( x*gamcoeff/(g_0*bar_xi) )^r;
            end
            %
            psi_n = double(limit( diff(gamma(2+s)^(tt-tau+1),s,n),s,-1 ));
            %
            sum_n = sum_n + psi_n/factorial(n) * sum_r;
        end
        sum_tau = sum_tau + prod_j/bar_xi * sum_n;
    end
    B = (B + sum_tau) * x*gamcoeff/g_0;
    %
    out = 1 - (A+B);
end