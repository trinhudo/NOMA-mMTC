function OP = OPS_PQoM_MTCD_II(M,PA_QoMS,CPA_QoMS,P_BPEH_t,varphi_t,noiseVar,R_M_QoMS,R_t_QoMS,tauBPEH,isQoMS)
    SINR1 = cell2mat( arrayfun( @(x) PA_QoMS(2)*P_BPEH_t(x,:).*varphi_t(x,:)...
                ./(CPA_QoMS(1)*P_BPEH_t(x,:).*varphi_t(x,:)+noiseVar),(1:(M-1))','UniformOutput',false) );
    SINR2 = cell2mat( arrayfun( @(x) CPA_QoMS(1)*P_BPEH_t(x,:).*varphi_t(x,:)./(noiseVar),(1:(M-1))','UniformOutput',false) );
    
    R1 = tauBPEH*log2(1+SINR1);
    R2 = tauBPEH*log2(1+SINR2);
    
    OP = 1-mean((R1 > R_M_QoMS).*(R2 > R_t_QoMS).*(isQoMS),2)./mean(isQoMS,2);
end