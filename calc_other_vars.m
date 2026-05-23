function [aa,mm,ss,ee,hh,p0,t0] = calc_other_vars(pp,rr,tt,uu)

global gamma ga gb gc gd ge gf gg gh gi gj % was USE PIGAMMA_PAR

% Everything is normalized, i.e. a_norm = a/u_ref = sqrt{gamma T_norm}
aa = gf*sqrt(tt);
mm = uu./aa;
ss = log(pp)-gamma*log(rr);
ee = gb*pp+.5d0.*rr.*uu.^2;
hh = ga*tt;
p0 = pp.*(1.+gd.*mm.^2).^ga;
t0 = tt.*(1.+gd.*mm.^2);

end