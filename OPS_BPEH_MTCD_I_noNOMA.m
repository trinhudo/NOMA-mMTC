function OP = OPS_BPEH_MTCD_I_noNOMA(M,P_BPEH_t,beta_t,phi_t,noiseVar,R_M_PQoMS,isEH)
    SINR = cell2mat( arrayfun( @(x) (1-beta_t(x+1).*isEH(x+1,:)).*P_BPEH_t(x,:).*phi_t(x+1,:)/noiseVar,(1:(M-1))','UniformOutput',false ) );
    SINR = [zeros(1,size(phi_t,2)) ; SINR];
    
    R = log2(1+SINR)/(M-1);
	OP= mean( R < R_M_PQoMS,2 ); 
end