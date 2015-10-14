% Solve shallow water equation
% nx - number of spatial grids in x
% ny - number of spatial grids in y
% cfl - CFL number
% t_F - final time
% hl - water height on the left of dam
% hr - water height on the right of dam
% rk -  1, Kinetic flux  2, Roe flux 
% fre - output frequency
% called by MATLAB command line:  kfvs_sw_dam(40, 40, 0.4, 4., 10., 5., 1, 100)
% Author: Luo Li, lluoac@ust.hk
function kfvs_sw_dam(nx,ny,cfl,t_F,hl,hr,rk,fre)
grav=9.8;
itmax=1000; % maximum of time steps
m=1; % number of ghost points in x
n=1; % number of ghost points in y
xmin=0.;
xmax=200.;
ymin=0.;
ymax=200.;
dx=(xmax-xmin)/nx; % space step length
dy=(ymax-ymin)/ny; % space step length
x=zeros(nx);
y=zeros(ny);
u0=zeros(nx+2*m,ny+2*n,3);  
u1=zeros(nx+2*m,ny+2*n,3);

for j=1:ny
    y(j)=ymin+0.5*dy+(j-1)*dy;
end

for i=1:nx
    x(i)=xmin+0.5*dx+(i-1)*dx;
end

% initialize
for j=1:ny
    for i=1:nx
        if x(i)<=0.5*(xmax-xmin)
           u0(m+i,n+j,1)=hl;
           u0(m+i,n+j,2)=0.0;
           u0(m+i,n+j,3)=0.0;		   
        else
           u0(m+i,n+j,1)=hr;
           u0(m+i,n+j,2)=0.0;
           u0(m+i,n+j,3)=0.0;	
        end
    end
end

% apply boundary condition
for j=1:2*n+ny
	for i=1:m
		% left outflow
		u0(i,j,:)=u0(m+i,j,:);
		% right outflow
		u0(m+nx+i,j,:)=u0(m+nx+i-m,j,:);
	end
end

for i=1:2*m+nx
	for j=1:n
		% bottom reflection
		u0(i,j,1)=u0(i,2*n+1-j,1);
		u0(i,j,2)=u0(i,2*n+1-j,2);
		u0(i,j,3)=-u0(i,2*n+1-j,3);
		% top reflection
		u0(i,n+ny+j,1)=u0(i,n+ny+1-j,1);
		u0(i,n+ny+j,2)=u0(i,n+ny+1-j,2);
		u0(i,n+ny+j,3)=-u0(i,n+ny+1-j,3);
	end
end

% dam
for j=1:ny
	if y(j)>=95. && y(j)<=170.
		u0(m+nx/2,n+j,:)  =u0(m+nx/2+2,n+j,:);
		u0(m+nx/2+1,n+j,:)=u0(m+nx/2-1,n+j,:);
	else 
		u0(m+nx/2,n+j,:)  =u0(m+nx/2-1,n+j,:);
		u0(m+nx/2+1,n+j,:)=u0(m+nx/2+2,n+j,:);
		u0(m+nx/2,n+j,2)  =-u0(m+nx/2-1,n+j,2);
		u0(m+nx/2+1,n+j,2)=-u0(m+nx/2+2,n+j,2);		
	end
end

surface(u0(m+1:m+nx,n+1:n+ny,1));
pause;

% start time stepping
t=0.0;
for it=1:itmax
	rhomaxx=0.;
	rhomaxy=0.;
	for j=1:ny
		for i=1:nx
			   % primitive variables
			   ux=u0(m+i,n+j,2)/u0(m+i,n+j,1);
			   uy=u0(m+i,n+j,3)/u0(m+i,n+j,1);
			   % sound speed
			   c=sqrt(grav*u0(m+i,n+j,1));
			   % evaluate the maximum of velocity
			   rhox = max(abs(ux+c),abs(ux-c));
			   rhomaxx = max(rhomaxx,rhox);
			   rhoy = max(abs(uy+c),abs(uy-c));
			   rhomaxy = max(rhomaxy,rhoy);			   
		end
	end
    % compute time step length 
    dt=cfl/max(rhomaxx/dx,rhomaxy/dy);
	
    % final time step length
    if t<t_F && t+dt>t_F
        dt=t_F-t;
    end
    t=t+dt
	
    % if reach time limit, stop time stepping
    if t+dt>t_F || it==itmax
 		surface(u0(m+1:m+nx,n+1:n+ny,1));       
        break;
    end
    
    lambdax=dt/dx;
    lambday=dt/dy;	
	
    % update solution
	if rk == 1	
	for j=1:ny
		for i=1:nx
			u1(m+i,n+j,:)=u0(m+i,n+j,:)+...
				lambdax*(k_f_sw(u0(m+i-1,n+j,:),u0(m+i,n+j,:))-k_f_sw(u0(m+i,n+j,:),u0(m+i+1,n+j,:)))+...
				lambday*(k_g_sw(u0(m+i,n+j-1,:),u0(m+i,n+j,:))-k_g_sw(u0(m+i,n+j,:),u0(m+i,n+j+1,:)));
		end
	end
	else
	for j=1:ny
		for i=1:nx
			u1(m+i,n+j,:)=u0(m+i,n+j,:)+...
				lambdax*(r_f_sw(u0(m+i-1,n+j,:),u0(m+i,n+j,:))-r_f_sw(u0(m+i,n+j,:),u0(m+i+1,n+j,:)))+...
				lambday*(r_g_sw(u0(m+i,n+j-1,:),u0(m+i,n+j,:))-r_g_sw(u0(m+i,n+j,:),u0(m+i,n+j+1,:)));
		end
	end
	end
	
	% apply boundary condition
	for j=1:2*n+ny
		for i=1:m
			% left outflow
			u1(i,j,:)=u1(m+i,j,:);
			% right outflow
			u1(m+nx+i,j,:)=u1(m+nx+i-m,j,:);
		end
	end

	for i=1:2*m+nx
		for j=1:n
			% bottom reflection
			u1(i,j,1)=u1(i,2*n+1-j,1);
			u1(i,j,2)=u1(i,2*n+1-j,2);
			u1(i,j,3)=-u1(i,2*n+1-j,3);
			% top reflection
			u1(i,n+ny+j,1)=u1(i,n+ny+1-j,1);
			u1(i,n+ny+j,2)=u1(i,n+ny+1-j,2);
			u1(i,n+ny+j,3)=-u1(i,n+ny+1-j,3);
		end
	end

	% dam
	for j=1:ny
		if y(j)>=95. && y(j)<=170.
			u1(m+nx/2,n+j,:)  =u1(m+nx/2+2,n+j,:);
			u1(m+nx/2+1,n+j,:)=u1(m+nx/2-1,n+j,:);
		else 
			u1(m+nx/2,n+j,:)  =u1(m+nx/2-1,n+j,:);
			u1(m+nx/2+1,n+j,:)=u1(m+nx/2+2,n+j,:);
			u1(m+nx/2,n+j,2)  =-u1(m+nx/2-1,n+j,2);
			u1(m+nx/2+1,n+j,2)=-u1(m+nx/2+2,n+j,2);		
		end
	end	  
	
    % swap
    u0(:,:,:)=u1(:,:,:);
	
    if mod(it,fre)==0 % draw frequency
		surface(u0(m+1:m+nx,n+1:n+ny,1));
		pause;
    end
	
end

fid = fopen('sw.plt', 'w');
fprintf(fid, 'TITLE     = \"example\"\n');
fprintf(fid, 'VARIABLES = \"x\"\t\"y\"\t\"h\"\t\"u\"\t\"v\"\n');
fprintf(fid, 'ZONE T=\"TTT\"\nI=%d, J=%d, ZONETYPE=Ordered\n',nx,ny);
fprintf(fid, 'DATAPACKING=POINT\n');
for j=1:ny
    for i=1:nx
        h=u0(m+i,n+j,1);
        ux=u0(m+i,n+j,2)/h;
        uy=u0(m+i,n+j,3)/h;		
		fprintf(fid, '%f \t%f \t%f \t%f \t%f\n',x(i),y(j),h,ux,uy);		
    end
end
fclose(fid);
