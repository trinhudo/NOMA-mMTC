% =================================================
%     m,n /   | (an,An) (ap,Ap) \
%    H    | z |                 | 
%     p,q \   | (bm,Bm) (bq,Bq) /
% =================================================
% =================================================
function out = foxH(an,An,ap,Ap...
                   ,bm,Bm,bq,Bq,z)
    %%
    warning('off','MATLAB:integral:MaxIntervalCountReached');
    %% Integrand definition
    F = @(s) (GammaProd(bm,Bm,s)...
           .* GammaProd(1-an,-An,s).*z.^-s)...
            ./ (GammaProd(1-bq,-Bq,s)...
                .* GammaProd(ap,Ap,s));
    %% Parameters:
    p = length([An Ap]);
    q = length([An Ap]);
    alphaFox = sum(An)-sum(Ap)+sum(Bm)-sum(Bq);
    mu = sum([Bm Bq]) - sum([An Ap]);
    betaFox = prod([An Ap].^-[An Ap])...
            * prod([Bm Bq].^-[Bm Bq]);
    delta = sum([bm bq]) - sum([an ap]) + (p-q)/2;
    %% Conditions per contour:
    % Contour L_(c+i*infinity):
    condition01= (alphaFox>0) && ...
        (abs(angle(z))<pi*alphaFox/2);
    condition02= (alphaFox==0)&& (delta*mu+...
        real(delta)) <-1 && (angle(z)==0);
    condition0 = condition01 || condition02;
    % contour L_(-infinity)
    condition11 = (mu>0) && (z~=0);
    condition12 = (mu==0) && (abs(z)<betaFox) && (abs(z)>0);
    condition13 = (mu==0) && (abs(z)==betaFox)...
        && (real(delta)<-1);
    condition1 = condition11||condition12||condition13;
    % contour L_(+infinity)
    condition21 = (mu<0) && (z~=0);
    condition22 = (mu==0) && (abs(z)>betaFox);
    condition2 = condition21 || condition22;
    %% Contour preparation:
    epsilon = 10^1.2;
    Sups = min((1-an)./An); Infs = max(-bm./Bm);
    if(isempty(Sups) && isempty(Infs))
        WPx=1;
    elseif(isempty(Sups) && ~isempty(Infs))
        WPx = Infs +epsilon;
    elseif(~isempty(Sups) && isempty(Infs))
        WPx = Sups -epsilon;
    else
        WPx = (Sups + Infs)/2;
    end
    WayPoints = [WPx-1i*epsilon WPx+1i*epsilon];
    %% integration:
    if(condition0 || (~condition1 && ~condition2))
        infity = 100;
        out = (1/(2i*pi))*integral(F,WPx-1i*infity,...
            WPx+1i*infity);
        return
    end
    if(~condition1)
        infity = 100;
        if(~isempty(min(-bm./Bm)))
            infity = infity - min(-bm./Bm);
        end
        out = (1/(2i*pi))*integral(F,-infity,...
            -infity,'Waypoints',WayPoints);
        return
    end
    if(condition2)
        infity = 100;
        if(~isempty(max((1-an)./An)))
            infity = infity + max((1-an)./An);
        end
        Tol = 10^-5;
        out = (1/(2i*pi))*integral(F,infity...
            ,infity,'Waypoints',WayPoints);
    end
%% ***** GammaProd subfunction *****
    function output = GammaProd(p,x,X)
        [pp,XX] = meshgrid(p,X);
        xx = meshgrid(x,X);
        if (isempty(p))
            output = ones(size(X));
        else
            output = reshape(prod(double(...
                gammaz(pp+xx.*XX)),2),size(X));
        end
    end
end