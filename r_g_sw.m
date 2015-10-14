% Roe flux for shallow water
% ul - 1x3 array, state on the bottom
% ur - 1x3 array, state on the top
% flux - 1x1x3 array, output flux
% Author Luo Li, lluoac@ust.hk
function [flux]=r_g_sw(ul,ur)
grav=9.8;
flux=zeros(1,1,3);
uavg=zeros(3);
lambda=zeros(3);
abslam=zeros(3);
diff=zeros(3);
diffusion=zeros(3);

% Roe averages
hl=ul(1);
hr=ur(1);
uxl=ul(2)/ul(1);
uxr=ur(2)/ur(1);
uyl=ul(3)/ul(1);
uyr=ur(3)/ur(1);

hlSqrt=sqrt(hl);
hrSqrt=sqrt(hr);

uavg(1)=.5*(hl+hr);
uavg(2)=(uxl*hlSqrt+uxr*hrSqrt)/(hlSqrt+hrSqrt);
uavg(3)=(uyl*hlSqrt+uyr*hrSqrt)/(hlSqrt+hrSqrt);

% sound speed
cavg=sqrt(grav*uavg(1));

lambda(1) = uavg(2)-cavg;
lambda(2) = uavg(2);
lambda(3) = uavg(2)+cavg;
abslam(1)=abs(lambda(1));
abslam(2)=abs(lambda(2));
abslam(3)=abs(lambda(3));

diff(1)=ur(1)-ul(1);
diff(2)=ur(2)-ul(2);
diff(3)=ur(3)-ul(3);

% diffusion
diffusion(1)=((abslam(1)*lambda(3) - lambda(1)*abslam(3))*diff(1)... 
			+ (abslam(3) - abslam(1))*diff(3))/(2.*cavg);
diffusion(2)=uavg(2)*((abslam(1)*lambda(3)...
			- lambda(1)*abslam(3))/(2.*cavg) - abslam(2))*diff(1)...
			+ abslam(2)*diff(2)+ uavg(2)*(abslam(3)-abslam(1))*diff(3)/(2.*cavg);
diffusion(3)=(lambda(1)*lambda(3)*(abslam(1)-abslam(3))*diff(1)... 
			+ (lambda(3)*abslam(3)-lambda(1)*abslam(1))*diff(3))/(2.*cavg);

% flux
flux(1,1,1) = .5*(uyl*hl+uyr*hr) - .5*diffusion(1); 	
flux(1,1,2) = .5*(uxl*uyl*hl+uxr*uyr*hr) - .5*diffusion(2);
flux(1,1,3) = .5*((uyl*uyl*hl+grav*hl*hl/2.)+(uyr*uyr*hr+grav*hr*hr/2.))... 
		- .5*diffusion(3);