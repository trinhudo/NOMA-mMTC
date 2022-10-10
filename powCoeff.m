function [PA,CPA] = powCoeff(K)
    % Power Allocation Array
    PA = (2.^(1:(K+1))-1)/(sum(2.^(1:(K+1))-1));
    % Cummulative Power Allocation, i.e., Portion of Intra-MTCD Interference Power
    CPA= [arrayfun(@(x) sum(PA(1:x)),K:-1:1),0];
    CPA(2:end) = fliplr(CPA(2:end));
end