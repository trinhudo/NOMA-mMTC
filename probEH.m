function out = probEH(M,Omega_t,PL_I2I)
    % Initialize
    rho_t = zeros(1,M);
    
    barXi = @(tau,t) prod(Omega_t(tau:t).*PL_I2I(tau:t));
    cF_Xi = @(xi,tau,t) meijerGtol([],[],[ones(1,t-tau),0],[],xi/barXi(tau,t));
    % Initial Value
    rho_t(1) = 0;
    rho_t(M) = 0;
    % Other Values
    for tt = 1:(M-2)
        rho_t(tt+1) = (1-rho_t(tt))*cF_Xi(1,tt+1,tt+1);
        for tau = 1:(tt-1)
            rho_t(tt+1) = rho_t(tt+1) + (1-rho_t(tau))*prod(rho_t((tau+1):tt)) * cF_Xi(1,tau+1,tt);
        end
    end
    % Output
    out = rho_t;
end