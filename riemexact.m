function [p,u,rhoc,rhod,el,er] = riemexact(aa,ab,pa,pb,ua,ub)

global gamma ga gb gc gd ge gf gg gh gi gj % era USE PIgamma_PAR

% Initial guess, left and right wave are expansion fans -> pseudo expansion
rhoa = gamma*pa/aa^2;
capaa = 1.0/ge/rhoa;
rhob = gamma*pb/ab^2;
capab = 1.0/ge/rhob;
elold = false;
erold = false;
el    = true;
er    = true;
pold = 0.5*(pa+pb);

while ((el ~= elold) || (er ~= erold))
    elold = el;
    erold = er;
    it = 0;
    f =1.0;
    pnew = pold *2.0;
    while (abs(f) > 1.e-15) && (abs((pnew-pold)/pnew) > 1.e-15)
        % treating it as a compression until I find the intersection
        if (it > 0)
            pold = pnew;
        end
        it = it+1;
        pp = pold + pold*1.e-5;
        pm = pold - pold*1.e-5;
        % Computes the difference between the left and right expansions
        % Newton Rhapson
        % Computes derivatives numerically

        % Lists all different cases, left expansion, right compressione etc
        % Also when the shock disappears
        if (el && er)
            fp = fel(pp,ua,aa,pa)      - fer(pp,ub,ab,pb);
            fm = fel(pm,ua,aa,pa)      - fer(pm,ub,ab,pb);
            f = fel(pold,ua,aa,pa)     - fer(pold,ub,ab,pb);
        elseif (el)
            fp = fel(pp,ua,aa,pa)      - fsr(pp,ub,capab,pb);
            fm = fel(pm,ua,aa,pa)      - fsr(pm,ub,capab,pb);
            f  = fel(pold,ua,aa,pa)    - fsr(pold,ub,capab,pb);
        elseif (er)
            fp = fsl(pp,ua,aa,pa)      - fer(pp,ub,ab,pb);
            fm = fsl(pm,ua,capaa,pa)   - fer(pm,ub,ab,pb);
            f  = fsl(pold,ua,capaa,pa) - fer(pold,ub,ab,pb);
        else
            fp = fsl(pp,ua,capaa,pa)   - fsr(pp,ub,capab,pb);
            fm = fsl(pm,ua,capaa,pa)   - fsr(pm,ub,capab,pb);
            f  = fsl(pold,ua,capaa,pa) - fsr(pold,ub,capab,pb);
        end
        
        ff = (fp-fm)/(pp-pm);
        
        pnew = pold - f/ff;
         
        if (pnew < 0.0d0)
            pnew = 1.e-6;
        end
    end
    
    p = pnew;
    if (el)
        u = ua + gg*aa*(1.d0-(p/pa)^gi);
    elseif (~el)
        u = ua - (p-pa)*sqrt(capaa/(p+pa/gc));
    end
    pc = p;
    pd = p;
    uc = u;
    ud = u;
    
    if (el)
        rhoc = rhoa*(pc/pa)^(1.0/gamma);
    elseif (~el)
        rhoc = rhoa*(1.0/gc+pc/pa)/(pc/pa/gc+1.0);
    end
    if (er)
        rhod = rhob*(pd/pb)^(1.0/gamma);
    elseif (~er)
        rhod = rhob*(1.0/gc+pd/pb)/(pd/pb/gc+1.0);
    end
    ac = sqrt(gamma*pc/rhoc);
    ad = sqrt(gamma*pd/rhod);
    
    ala = ua-aa;
    alc = uc-ac;
    
    ald = ud+ad;
    alb = ub+ab;
    
    el = false;
    er = false;
    % Checks wheter the waves converge  and computes everything
    if (ala <= alc) 
        el = true;
    end
    if (ald <= alb)
        er = true;
    end
    
    pold = pc;
end

    function fel = fel(p,ua,aa,pa)
        %      use pigamma_par
        fel = ua - gg*aa*((p/pa)^gi-1.0);
    end

    function fer = fer(p,ub,ab,pb)
        %      use pigamma_par
        fer = ub + gg*ab*((p/pb)^gi-1.0);
    end

    function fsl = fsl(p,ua,capaa,pa)
        %      use pigamma_par
        fsl = ua - (p-pa)*sqrt(capaa/(p+pa/gc));
    end

    function fsr = fsr(p,ub,capab,pb)
        %      use pigamma_par
        fsr = ub + (p-pb)*sqrt(capab/(p+pb/gc));
    end

end
% ************************************************************************      