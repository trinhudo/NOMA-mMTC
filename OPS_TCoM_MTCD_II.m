function OP = OPS_TCoM_MTCD_II(M,K_t,PA_CoMS,CPA_CoMS,P_BTEH_t,varphi_tk,noiseVar,R_M_CoMS,R_tk_CoMS,tauBTEH,isCoMS)

    OP = ones(M-1,max(K_t));
    for tt = 1:M-1
        
%         a_1 = 2^( (M-1)*R_M_CoMS )-1;
%         a_2 = 2^( (M-1)*R_M_CoMS/(1-alpha_t(tt+1)) )-1;
        
        for kk = 1:K_t(tt)
            varphi = squeeze(varphi_tk(tt,kk,:))';
            isExist= squeeze(isCoMS(tt,kk,:))';
            
            sinr18 = PA_CoMS(tt,K_t(tt)+1)*P_BTEH_t(tt,:).*varphi...
                ./(CPA_CoMS(tt,1)*P_BTEH_t(tt,:).*varphi+noiseVar);
            rate18 = tauBTEH(tt+1,:).*log2( 1+sinr18 );
            
            SPEvent = (rate18 >= R_M_CoMS);
            
%             qqs = (kk:K_t(tt))+1;
%             b_1 = 2^( (M-1)*R_tk_CoMS(tt,kk) )-1;
%             b_3 = max( [ a_1/(PA_CoMS(tt,K_t(tt)+1) - a_1*CPA_CoMS(tt,1)),...
%                         b_1./(PA_CoMS(tt,kk:K_t(tt))- b_1*CPA_CoMS(tt,qqs))] );                   
%             b_2 = 2^( (M-1)*R_tk_CoMS(tt,kk)/(1-alpha_t(tt+1)) )-1;
%             b_4 = max( [ a_2/(PA_CoMS(tt,K_t(tt)+1) - a_2*CPA_CoMS(tt,1)),...
%                         b_2./(PA_CoMS(tt,kk:K_t(tt))- b_2*CPA_CoMS(tt,qqs))] );
%             
%             SPEvent21= (P_BTEH_t(tt,:).*varphi/noiseVar >= b_3).*(isEH(tt+1,:)==0);
%             SPEvent22= (P_BTEH_t(tt,:).*varphi/noiseVar >= b_4).*(isEH(tt+1,:)==1);
%             OP(tt,kk) = 1 - mean(SPEvent.*SPEvent21.*isExist)/mean(isExist)...
%                           - mean(SPEvent.*SPEvent22.*isExist)/mean(isExist);
            for qq = kk:K_t(tt)
                sinr19 = PA_CoMS(tt,qq)*P_BTEH_t(tt,:).*varphi...
                    ./(CPA_CoMS(tt,qq+1)*P_BTEH_t(tt,:).*varphi+noiseVar);
                rate19 = tauBTEH(tt+1,:).*log2( 1+sinr19 );
                
                SPEvent = SPEvent .* (rate19 >= R_tk_CoMS(tt,qq));
            end
            OP(tt,kk) = 1 - mean(SPEvent.*isExist)/mean(isExist);
        end
    end
end