function out = OPe2e_XQoM_MTCD_II(OP_XQoM_MTCD_II,M,OP_BXEH_MTCD_I)
    SPe2e_XQoM_MTCD_II = ones(size(OP_XQoM_MTCD_II));
    for tt = 1:(M-1)
        %_
        SP_BXEH_MTCD_I = ones(1,size(OP_BXEH_MTCD_I,2));
        if (tt >= 2)
            for ii = 1:tt
                SP_BXEH_MTCD_I = SP_BXEH_MTCD_I .* (1-OP_BXEH_MTCD_I(tt+1,:));
            end
        end
        %_
        SPe2e_XQoM_MTCD_II(tt,:) = SP_BXEH_MTCD_I.*(1-OP_XQoM_MTCD_II(tt,:));
    end
    out = 1 - SPe2e_XQoM_MTCD_II;
end