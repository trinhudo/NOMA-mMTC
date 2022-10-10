function OP = OPS_BTEH_MTCD_I(M,PA_QoMS,CPA_QoMS,P_BTEH_t,phi_t,noiseVar,R_M_TQoMS,isExist,tauBTEH)
	SINR = cell2mat( arrayfun( @(x) PA_QoMS(2)*P_BTEH_t(x,:).*phi_t(x+1,:)...
                ./(CPA_QoMS(1)*P_BTEH_t(x,:).*phi_t(x+1,:)+noiseVar),(1:(M-1))','UniformOutput',false) ).*isExist...
         + cell2mat( arrayfun( @(x) P_BTEH_t(x,:).*phi_t(x+1,:)/noiseVar,(1:(M-1))','UniformOutput',false) ).*(~isExist);
    SINR = [zeros(1,size(phi_t,2)); SINR];
    
	R = tauBTEH.*log2(1+SINR);
	OP= mean( R < R_M_TQoMS,2 ); 
end