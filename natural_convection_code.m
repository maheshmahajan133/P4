%% finding the temp. variation for the solid part
 
% given parameters
dx=0.02;dy=dx;
D=1.0;
nx=(D/dx)+1;
ny=(D/dy)+1;
mu=1.8e-5;cp=1005;
Tc=20;q=100;k=mu*cp/0.72;h=2;
iter=1000;
Tw=(q/h)+Tc;
 
T=zeros(nx,ny);
 
%boundary conditions for the solid part
T(nx,:)=Tw;
T(41,:)=T(42,:)+q*dx/k;
T(:,1)=T(:,2);
T(:,ny)=T(:,ny-1);
 
 
for n=1:iter
for i=41:nx-1
    for j=2:ny-1
        
        T(i,j)=0.25*(T(i-1,j)+T(i+1,j)+T(i,j-1)+T(i,j+1));
        T(41,:)=T(42,:)+q*dx/k;
        T(:,1)=T(:,2);
        T(:,ny)=T(:,ny-1);
    end
end
end
 
% normalizing the temp.
for i=41:nx
    for j=1:ny
        T(i,j)=(T(i,j)-Tc)*k/(D*q);
    end
end
 
% finding the temp. at the cell midpoint(for further calculations in fluid
% part)
for j=2:ny
    T(41,j)=(T(41,j)+T(41,j-1))/2;
end
 
 
 
 
%% finding the temp. variation of fluid part
nx=41; ny=51; D=1.0;
dx=(D)/(nx-1); dy=D/(ny-1); dt=0.00001; %Re=10000;
un=0.0; us=0.0; gamma=1.0; Pr=0.72; kr=1.0; Ra=100000;
maxit=20000; phi=30;
 
% '_tilda' denotes assumed value
% '_n' denotes corrected value
u_tilda=0.0*ones(nx,ny+1); u_n=0.0*ones(nx,ny+1);
v_tilda=zeros(nx+1,ny); v_n=zeros(nx+1,ny);
p_dash=0.0*ones(nx+1,ny+1); p_n=0.0*ones(nx+1,ny+1);
T_n_plus_1=0.0*ones(nx+1,ny+1); % Temp. at timestep n
T_n=0.0*ones(nx+1,ny+1); % temp. at timestep n+1
 
% values at normal grid point
u_cell=zeros(nx,ny); v_cell=zeros(nx,ny); p_cell=zeros(nx,ny); T_cell=zeros(nx,ny);
 
nx=41;
for t=0:maxit
    %B.Cs of u velocity
    for i=1:nx
        u_n(i,1)=(2*us)-u_n(i,2);
        u_n(i,ny+1)=(2*un)-u_n(i,ny);
    end
    for j=1:ny+1
        u_n(1,j)=0.0;
        u_n(nx,j)=0.0;
    end
    
    %B.Cs of v velocity
    for j=1:ny
        v_n(1,j)=-v_n(2,j);
        v_n(nx+1,j)=-v_n(nx,j);
    end
    for i=1:nx+1
        v_n(i,1)=0.0;
        v_n(i,ny)=0.0;
    end
    
    %solving the temperature transport equation
    for i=2:nx
        for j=2:ny
            axd=dx/dy; ayd=dy/dx;
            aec=u_n(i,j)*dy;
            awc=u_n(i-1,j)*dy;
            anc=v_n(i,j)*dx;
            asc=v_n(i,j-1)*dx;
            
            %hybrid scheme
            ae=max([-aec,axd-aec/2,0]);
            aw=max([awc,axd+awc/2,0]);
            an=max([-anc,axd-anc/2,0]);
            as=max([asc,axd+asc/2,0]);
            ap=ae+aw+an+as+(aec-awc+anc-asc);
            
           
            T_n_plus_1(i,j)=T_n(i,j)...
                -dt/dx/dy*(ap*T_n(i,j)-ae*T_n(i+1,j)-aw*T_n(i-1,j)-an*T_n(i,j+1)-as*T_n(i,j-1));
        end
    end
    
    %B.Cs of temperature
    for j=1:ny
        T_n_plus_1(1,j)=(-1*dx)+T_n_plus_1(2,j);
        T_n_plus_1(nx+1,j)=2*T(41,j)-T_n_plus_1(nx,j);
    end
    
    for i=1:nx+1
        T_n_plus_1(i,1)=T_n_plus_1(i,2);
        T_n_plus_1(i,ny+1)=T_n_plus_1(i,ny);
    end
    
    T_n=T_n_plus_1;
    
    % solving for x momentum equation
    for i=2:nx-1
        for j=2:ny
            
            axd=gamma*Pr*dy/dx; ayd=gamma*Pr*dx/dy;  
            aec=u_n(i+1,j)*dy;
            awc=u_n(i-1,j)*dy;
            anc=(v_n(i,j)+v_n(i+1,j))*0.5*dx;
            asc=(v_n(i,j-1)+v_n(i+1,j-1))*0.5*dx;
            
            %hybrid scheme
            ae=max([-aec,axd-aec/2,0]);
            aw=max([awc,axd+awc/2,0]);
            an=max([-anc,axd-anc/2,0]);
            as=max([asc,axd+asc/2,0]);
            ap=ae+aw+an+as+(aec-awc+anc-asc);
            
                u_tilda(i,j)=u_n(i,j)...
                    -(dt/dx)*(p_n(i+1,j)-p_n(i,j))...
                    -(dt/dx/dy)*(ap*u_n(i,j)-ae*u_n(i+1,j)-aw*u_n(i-1,j)-an*u_n(i,j+1)-as*u_n(i,j-1))...
                    +Ra*Pr*cosd(phi)*0.5*(T_n_plus_1(i+1,j)+T_n_plus_1(i,j))*dt;
        end
    end
    for i=1:nx
        u_tilda(i,1)=(2*us)-u_tilda(i,2);
        u_tilda(i,ny+1)=(2*un)-u_tilda(i,ny);
    end
    for j=1:ny+1
        u_tilda(1,j)=0.0;
        u_tilda(nx,j)=0.0;
    end
    
    % solving for y momentum equation
    for i=2:nx
        for j=2:ny-1
            axd=gamma*Pr*dy/dx; ayd=gamma*Pr*dx/dy; 
            aec=(u_n(i,j)+u_n(i,j+1))*0.5*dy;
            awc=(u_n(i-1,j)+u_n(i-1,j+1))*0.5*dy;
            anc=v_n(i,j+1)*dx;
            asc=v_n(i,j-1)*dx;
            
            %hybrid scheme
            ae=max([-aec,axd-aec/2,0]);
            aw=max([awc,axd+awc/2,0]);
            an=max([-anc,axd-anc/2,0]);
            as=max([asc,axd+asc/2,0]);
            ap=ae+aw+an+as+(aec-awc+anc-asc);
                v_tilda(i,j)=v_n(i,j)...
                    -(dt/dx)*(p_n(i,j+1)-p_n(i,j))...
                     -(dt/dx/dy)*(ap*v_n(i,j)-ae*v_n(i+1,j)-aw*v_n(i-1,j)-an*v_n(i,j+1)-as*v_n(i,j-1))...
                    +Ra*Pr*sind(phi)*0.5*(T_n_plus_1(i,j+1)+T_n_plus_1(i,j))*dt;
        end
    end
    for j=1:ny
        v_tilda(1,j)=-v_tilda(2,j);
        v_tilda(nx+1,j)=-v_tilda(nx,j);
    end
    for i=1:nx+1
        v_tilda(i,1)=0.0;
        v_tilda(i,ny)=0.0;
    end
    
    % solving for pressure correction
    p_dash=p_n; error=10.0; max1=0.0; iter=0;
    while(error>0.0001)
        for i=2:nx
            for j=2:ny
                    p_dash(i,j)=(0.5/(dx^2+dy^2))*(((dy^2)*(p_dash(i+1,j)+p_dash(i-1,j)))...
                        +((dx^2)*(p_dash(i,j+1)+p_dash(i,j-1)))...
                        -((dx*dy/dt)*(((u_tilda(i,j)-u_tilda(i-1,j))*dy)...
                        +((v_tilda(i,j)-v_tilda(i,j-1))*dx))));
            end
        end
        max1=0.0;
        for i=2:nx
            for j=2:ny
                error=abs(p_dash(i,j)-p_n(i,j));
                if (error>max1)
                    max1=error;
                end
            end
        end
        error=max1;
        p_n=p_dash;
        iter=iter+1;
    end
    
    % finding the corrected u
    for i=2:nx-1
        for j=2:ny
            u_n(i,j)=u_tilda(i,j)-((dt/dx)*(p_dash(i+1,j)-p_dash(i,j)));
        end
    end
    
    % finding the corrected v
    for i=2:nx
        for j=2:ny-1
            v_n(i,j)=v_tilda(i,j)-((dt/dy)*(p_dash(i,j+1)-p_dash(i,j)));
        end
    end
    
end
 
 
% calculating the values at cell node points
for i=1:nx
    for j=1:ny
        u_cell(i,j)=(u_n(i,j+1)+u_n(i,j))/2;
        v_cell(i,j)=(v_n(i+1,j)+v_n(i,j))/2;
        p_cell(i,j)=(p_n(i,j)+p_n(i+1,j)+p_n(i+1,j+1)+p_n(i,j+1))/4;
        T_cell(i,j)=(T_n_plus_1(i,j)+T_n_plus_1(i+1,j)+T_n_plus_1(i+1,j+1)+T_n_plus_1(i,j+1))/4;
    end
end
 
% w=zeros(nx,ny);
% for i=2:nx-1
%     for j=2:ny-1
%         w(i,j)=(-(u_cell(i,j+1)-u_cell(i,j-1))/(2*dy))+((v_cell(i+1,j)-v_cell(i-1,j))/(2*dx));
%     end
% end
 
% appending the covection obtained temp. values to the global temp. matrix
for i=1:41
    for j=1:51
        T(i,j)=T_cell(i,j);
    end
end
 
T_trans=T';
 
for i=1:51
    T1(i,:)= T_trans(52-i,:);
end

% plotting the stramlines and isotherms for the u,v and T values obtained. 

[x,y]=meshgrid(linspace(0,1,5),linspace(0,1,5));
lines=streamline(stream2(linspace(0,1,nx),linspace(0,1,ny),u_cell',v_cell',x,y));
contour(T1)
