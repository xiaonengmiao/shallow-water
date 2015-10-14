% Kinetic flux for shallow water
% ul - 1x3 array, state on the bottom
% ur - 1x3 array, state on the top
% flux - 1x1x3 array, output flux
% Author Luo Li, lluoac@ust.hk
function [flux]=k_g_sw(ul,ur)
grav=9.8;
flux=zeros(1,1,3);
% bottom state
hl=ul(1);
ulx=ul(2)/ul(1);
uly=ul(3)/ul(1);
laml=1./grav/hl;
u0l=0.5*erfc(-sqrt(laml)*uly);
u1l=uly*u0l+0.5*exp(-laml*uly^2)/sqrt(pi*laml);
u2l=uly*u1l+0.5/laml*u0l;

% top state
hr=ur(1);
urx=ur(2)/ur(1);
ury=ur(3)/ur(1);
lamr=1./grav/hr;
u0r=0.5*erfc(sqrt(lamr)*ury);
u1r=ury*u0r-0.5*exp(-lamr*ury^2)/sqrt(pi*lamr);
u2r=ury*u1r+0.5/lamr*u0r;

% flux
flux(1,1,1)=hl*u1l*1.*1.+hr*u1r*1.*1.;
flux(1,1,2)=hl*u1l*ulx+hr*u1r*urx;
flux(1,1,3)=hl*u2l*1.*1.+hr*u2r*1.*1.;

