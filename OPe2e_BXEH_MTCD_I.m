function out = OPe2e_BXEH_MTCD_I(M,OP_BXEH_MTCD_I)
    SPe2e_BXEH_MTCD_I = ones(1,size(OP_BXEH_MTCD_I,2)); 
    for tt = 1:(M-1)
        SPe2e_BXEH_MTCD_I = SPe2e_BXEH_MTCD_I.*(1-OP_BXEH_MTCD_I(tt+1,:));
    end
    out = 1 - SPe2e_BXEH_MTCD_I;
end