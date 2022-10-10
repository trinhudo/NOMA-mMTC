% =================================================
%     m,n /   | (an,An) (ap,Ap) \
%    H    | z |                 |  z > 0
%     p,q \   | (bm,Bm) (bq,Bq) /
% =================================================
% =================================================
function out = WSCPfoxH(m,n,p,q,P,Q,z)
    %_
    if (n > 0)
        an = P(1:n,1)'; An = P(1:n,2)';
    else
        an = []; An = []; 
    end
    
    if (m > 0)        
        bm = Q(1:m,1)'; Bm = Q(1:m,2)';
    else
        bm = []; Bm = []; 
    end
    
    if (p-n > 0)    
        ap = P((n+1):p,1)'; Ap = P((n+1):p,2)';
    else  
        ap = []; Ap = [];
    end
    
    if (q-m > 0)
        bq = Q((m+1):q,1)'; Bq = Q((m+1):q,2)';
    else
        bq = []; Bq = [];
    end
    %_
	out = real(foxH(an,An,ap,Ap,bm,Bm,bq,Bq,z));
end