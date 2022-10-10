function OP = OPS_BTEH_MTCD_I_noNOMA(M,P_BTEH_t,phi_t,noiseVar,R_M_TQoMS,tauBTEH)
	SINR = cell2mat( arrayfun( @(x) P_BTEH_t(x,:).*phi_t(x+1,:)/noiseVar,(1:(M-1))','UniformOutput',false) );
    SINR = [zeros(1,size(phi_t,2)); SINR];
    
	R = tauBTEH.*log2(1+SINR);
	OP= mean( R < R_M_TQoMS,2 ); 
end