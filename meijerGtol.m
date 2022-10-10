function out1 = meijerGtol(An,Aq,Bm,Bp,x)
    if (x > 0)
        if isempty(An) && isempty(Aq) && isempty(Bp)
            if (Bm == 0)
                out1 = exp(-x);
            elseif (isequal(Bm,[1,0]))
                out1 = 2*x^(1/2)*besselk(1, 2*x^(1/2));
             elseif (isequal(Bm,[1,1,0]))
                if (x >= 2.1e3)
                    out1 = 0;
                else
                    out1 = meijerGHC(An,Aq,Bm,Bp,x);
                end
            else
                if (x >= 5e3)
                    out1 = 0;
                else
                    out1 = meijerGHC(An,Aq,Bm,Bp,x);
                end
            end
        else
            out1 = meijerGHC(An,Aq,Bm,Bp,x);
        end
    else
        out1 = 0;
    end
function out2 = meijerGHC(An,Ap,Bm,Bq,x)
    infty = 1e3;
    
    m = length(Bm);
    n = length(An);
    p = length(An) + length(Ap);
    q = length(Bm) + length(Bq);
    
    result = WSCPfoxH(m,n,p,q,[An', ones(n,1);Ap' ones(p-n,1)],...
        [ Bm' ones(m,1); Bq' ones(q-m,1) ], x);
%     meijerG(An,Ap,Bm,Bq,x);
    if (result < 0) || isnan(result) || (abs(result) >= infty)
        digits(64);  % fixed precision
        out2 = double(vpa(meijerG(An,Ap,Bm,Bq,sym(x))));
        digits(16);  % default precision
    else
        out2 = double(result);
    end