function OP = OPS_PCoM_MTCD_II(M,K_t,PA_CoMS,CPA_CoMS,P_BPEH_t,varphi_tk,noiseVar,R_M_CoMS,R_tk_CoMS,tauBPEH,isCoMS)

    OP = ones(M-1,max(K_t));
    for tt = 1:M-1
        for kk = 1:K_t(tt)
            varphi = squeeze(varphi_tk(tt,kk,:))';
            isExist= squeeze(isCoMS(tt,kk,:))';
            
            sinr18 = PA_CoMS(tt,K_t(tt)+1)*P_BPEH_t(tt,:).*varphi...
                ./(CPA_CoMS(tt,1)*P_BPEH_t(tt,:).*varphi+noiseVar);
            rate18 = tauBPEH*log2( 1+sinr18 );
            
            SPEvent = (rate18 >= R_M_CoMS);
            for qq = kk:K_t(tt)
                sinr19 = PA_CoMS(tt,qq)*P_BPEH_t(tt,:).*varphi...
                    ./(CPA_CoMS(tt,qq+1)*P_BPEH_t(tt,:).*varphi+noiseVar);
                rate19 = tauBPEH*log2( 1+sinr19 );
                
                SPEvent = SPEvent .* (rate19 >= R_tk_CoMS(tt,qq));
            end
            OP(tt,kk) = 1-mean(SPEvent.*isExist)/mean(isExist);
        end
    end
end