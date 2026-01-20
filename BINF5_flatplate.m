%FVM FLAT PLATE WITH INFLATION LAYER IN BOUNDARY. LEADING AND TRAILING EDGE
%FREE FLOW UPSTREAM AND DOWNSTREAM ZONES
%SIMULATED WITH SYMMETRIC  BOUNDARY  CONDITION
%ALGORITHM: SIMPLE
%SCHEME: UPWIND
% MATRIX SOL. METHOD: GAUSS
%FLUID: AIR AT 20° AND SEA LEVEL
%FACE INTERPOLATION METHOD: RIE-CHOW
%CO-LOCATED GRID WITHOUT GHOST CELLS
%STEADY STATE INCOMPRESSIBLE NAVIER STOKES

%clear all
close all

%% DATA

%Geometry
platLgt=0.6;                    %Flat plate length (m)
domHgt=1;                       %Domain Height (m)
domLgt=1.8;                     %Domain Length (m)
distPlat=0.5*(domLgt -platLgt); %Distance from origin to leading edge

%Fluid

rho=1.2;                %Density (Kg/m3)
mu=1.8e-5;              %Molecular Viscosity (N*s/m^2)
nu =mu/rho;             %Dinamic Viscosity (m/s^2)
u0=0.0125;              %Velocity at the inlet (m/s) /max reached=0.07; 0.0125;    
p0=1;                   %Outlet pressure (Prescribed)
Re=u0*rho*platLgt/mu;   %Reynolds number

%% GRID

%FLAT PLATE MESH 
%CODE  TO GENERATE MESH OVER A PLATE WITH INFLATION LAYERS AT THE 
%BOUNDARY, LEADING AND TRAILING EDGES 
%THERE IS A FREE - STREAM ZONE AT THE LEADING EDGE AND OTHER IN THE 
%TRAILING EDGE 

if (distPlat+ platLgt)<domLgt
    %VALID GEOMETRY
    
    %NUMBER OF CELLS IN EACH SECTION
    nx=93;% number of cells in x axis 
    ny=38;% number of cells in y axis
    ncells=nx*ny; %Total cells

    x=zeros(1,nx+1); %grid positions for x
    y=zeros(1,ny+1); %grid positions for y

    %Calculate the ratio between the length of the plate and the
    %  domain length
    fpRatio =platLgt/domLgt;

    %Calculate the ratio between the distance to the leading edge and the
    %  domain lenght origin (we will call this upstream zone)
    
    upstrRatio=distPlat/domLgt;


    %Divide proportionaly the number of cells in x between the plate 
    % and the remaining space 

    nx_fp=round(nx*fpRatio,0);%cell in the plate 
    nx_upstr=round(nx*upstrRatio);%Cells in the upstream zone
    nx_dwnstr=nx-(nx_fp + nx_upstr);%Cells in the downstream zone

    nxSymcond=[1:1:nx_upstr,nx-nx_dwnstr+1:1:nx];%indexes for symetric 
    %boundary condition

    %INFLATION LAYERS

    %INFLATION LAYER FOR X AXIS 

    dx_1=0.0033;%Lenght of the first cell at the leading and trailing edge
    grthFx=1.3;% Geometrical growth  Factor
    ninfLx=8;%Cells in the l/t edges side of the inflation layers

    %Calculate grid spaces for x in the inflation layer

    %define where  finishes the upstream zone and begins the plate
    x(nx_upstr +1)=distPlat;

    %Define where finishes the flat plate
    x(nx_upstr+nx_fp+1)=distPlat+ platLgt;

    %Define where finishes the domain
    x(nx+1)=domLgt;

    %Calculate grid spaces for x in the inflation layer (leading edge ,RS)
    %At the same time mirror the steps for the free stream portion
    %repeat for the trailing edge
    for i=1:ninfLx
        dx_inf=dx_1*grthFx^(i-1);%Layer lenght
        x(nx_upstr+1+i)=dx_inf + x(nx_upstr+i);%Flat plate LE
        x(nx_upstr+1-i)=-dx_inf +x(nx_upstr+2-i);%upwind portion
        
        %Apply same for the trailing edge

        x(nx_upstr + nx_fp + 1 + i) = dx_inf + ...
            x(nx_upstr + nx_fp + i);%Flat plate TE
        x(nx_upstr + nx_fp + 1-i) = -dx_inf + ...
            x(nx_upstr + nx_fp + 2-i);%Downstream portion
    end

    %SMOTH TRANSITIONS FOR X

    %SMOTH TRANSITION IN THE UPSTREAM ZONE

    dxMaxinf=x(nx_upstr + ninfLx +1)- x(nx_upstr + ninfLx); % Higher x grid space 
    % in inflation layers for X

    ninfLxupstr= nx_upstr - ninfLx;%Cells in the smooth transition 
    % upstream layer

    dxT2upstr=x(ninfLxupstr+1);%Lenght of the 
    %upstream smooth transition layer

    %%% Computig geometrical growth factor for upstream smooth 
    % transition layer

    grthFT2upstr=0.00001;%Geometrical grow factor 2
    dxT2upstrit=0;%Boundary layer height iterated
    while dxT2upstrit<=dxT2upstr
        grthFT2upstr=grthFT2upstr+0.00001;
        dxT2upstrit=dxMaxinf*((1-grthFT2upstr^ninfLxupstr)/...
            (1-grthFT2upstr));
    end

    %Fill the x grid vector with smoth transition downstream layer steps 
    
    for i=1:ninfLxupstr-1
        dx_inf=dxMaxinf*grthFT2upstr^(i-1);%layer lenght
        x(nx_upstr-ninfLx+1-i) = x(nx_upstr-ninfLx+2-i) - dx_inf;
    end


    %SMOTH TRANSITION IN THE DOWNSTREAM ZONE 

    ninfLxdwnstr=nx_dwnstr-ninfLx;%Cells in the smooth transition 
    % downstream layer

    dxT2dwnstr=domLgt-x(nx_upstr+nx_fp+ninfLx+1);%Lenght of the 
    %downstream smooth transition layer

    %%% Computig geometrical growth factor for downstream smooth 
    % transition layer

    grthFT2dwnstr=0.00001;%Geometrical grow factor 2
    dxT2dwnstrit=0;%Boundary layer height iterated
    while dxT2dwnstrit<=dxT2dwnstr
        grthFT2dwnstr=grthFT2dwnstr+0.00001;
        dxT2dwnstrit=dxMaxinf*((1-grthFT2dwnstr^ninfLxdwnstr)/...
            (1-grthFT2dwnstr));
    end

    %Fill the x grid vector with smoth transition downstream layer steps 
    
    for i=1:ninfLxdwnstr
        dx_inf=dxMaxinf*grthFT2dwnstr^(i-1);%layer lenght
        x(nx_upstr+nx_fp+ninfLx+1+i) = x(nx_upstr+nx_fp+ninfLx+i) + dx_inf;
    end


    %SMOTH TRANSITION IN THE PLATE

    %Calculate the distance from the origin to the midpoint of the plate
    distMidPlat=distPlat + 0.5*platLgt;

    %Smoth transition from leading edge inflation last layer to mid point
    %of the plate;

    %cells in the smooth transition zone of the plate (mid region)
    ninfLxmidFp=nx_fp-2*ninfLx;

    %Divide this distance in a upstream and downstream zone
    ninfLxmidFupstr= round(0.5*ninfLxmidFp);% Cells in the upstream zone

    ninfLxmidFdwnstr=ninfLxmidFp-ninfLxmidFupstr;%Cells in the downstream z.

    %SMOTH TRANSITION IN THE UPSTREAM ZONE (PLATE)

    dxT2upstrmidPlat=distMidPlat - x(nx_upstr + ninfLx + 1);%Lenght of the 
    %upstream smooth transition layer

    %%% Computig geometrical growth factor for upstream smooth 
    % transition layer

    grthFT2upstrmidPlat=0.00001;%Geometrical grow factor 2
    dxT2upstritmidPlat=0;%Boundary layer height iterated
    while dxT2upstritmidPlat<=dxT2upstrmidPlat
        grthFT2upstrmidPlat=grthFT2upstrmidPlat+0.00001;
        dxT2upstritmidPlat=dxMaxinf*((1-grthFT2upstrmidPlat^ninfLxmidFupstr)/...
            (1-grthFT2upstrmidPlat));
    end

    %Fill the x grid vector with smoth transition upstream layer steps 
    %for the plate 
    
    for i=1:ninfLxmidFupstr
        dx_inf=dxMaxinf*grthFT2upstrmidPlat^(i-1);%layer lenght
        x(nx_upstr + ninfLx + 1 + i) = x(nx_upstr + ninfLx + i)+ dx_inf;
    end

    %SMOTH TRANSITION IN THE DOWNSTREAM ZONE (PLATE)

    %Smoth transition from trailing edge inflation last layer to mid point
    %of the plate;

    %Lenght of the downstream smooth transition layer

    dxT2dwnstrmidPlat=x(nx - (nx_dwnstr + ninfLx) +1) - distMidPlat;

    %%% Computig geometrical growth factor for downstream smooth 
    % transition layer on plate

    grthFT2dwnstrmidPlat=0.00001;%Geometrical grow factor 2
    dxT2dwnstritmidPlat=0;%Boundary layer height iterated
    while dxT2dwnstritmidPlat<=dxT2dwnstrmidPlat
        grthFT2dwnstrmidPlat=grthFT2dwnstrmidPlat+0.00001;
        dxT2dwnstritmidPlat=dxMaxinf*((1-grthFT2dwnstrmidPlat^ninfLxmidFdwnstr)/...
            (1-grthFT2dwnstrmidPlat));
    end

    %Fill the x grid vector with smoth transition downstream layer steps 
    %for the plate 
    
    for i=1:ninfLxmidFdwnstr-1
        dx_inf=dxMaxinf*grthFT2dwnstrmidPlat^(i-1);%layer lenght
        x( nx - (nx_dwnstr + ninfLx) +1 - i) = x(nx - ...
            (nx_dwnstr + ninfLx) +2 - i)- dx_inf;
    end

    dxVec=diff(x);%X grid spaces vector

    %______________   Inflation layer parameters for Y axis  _____________

    dy_1=0.00167*0.25;% Height of the first cell (m);0.00021,0.0001;%
    grthF=1.15;% Geometrical growth  Factor ,1.3
    ninfL=10;% Cells in the boundary inflation layer
    
    
    %Calculate grid spaces for y in the inflation layer
    for i=1:ninfL
        y(i+1)=dy_1*grthF^(i-1) + y(i);
    end
    
    %Calculating inflation layer parameters to fill the remaining
    %spaces in  y grid vector and have a smooth transition
    
    dy_2=y(ninfL+1)-y(ninfL); % Higher y grid space in boundary inf layer
    ninfL2=ny-ninfL;%Cells in the second inflation layer
    dyt2=domHgt - y(ninfL+1);%Height of the second inflation layer
    
    %%% Computig geometrical growth factor for second inflation layer
    grthF2=0.00001;%Geometrical grow factor 2
    dyt2it=0;%Boundary layer height iterated
    while dyt2it<=dyt2
        grthF2=grthF2+0.00001;
        dyt2it=dy_2*((1-grthF2^ninfL2)/(1-grthF2));
    end
    
    
    % Fill the remaining y grid vector with second inflation layer
    for i=ninfL+1:ny
        y(i+1) = y(i) + dy_2*grthF2^(i-ninfL-1);
    end
    
    dyVec=diff(y);%Y grid spaces vector
    dyVec=flip(dyVec);
    
    y=flip(y);
    
    [X,Y]=meshgrid(x,y);
    
    %Calculate cell center coordinates
    
    xctrs=zeros(1,nx);% cell center x Grid vector
    yctrs=zeros(1,ny);% cell center y grid vector 
    
    for i=2:nx+1
        xctrs(i-1)=x(i)-0.5*dxVec(i-1);
    end
    
    for i=1:ny
        yctrs(i)=y(i)-0.5*dyVec(i);
    end
    
    [Xctrs,Yctrs]=meshgrid(xctrs,yctrs); %Cell centers coordinates

    %Calculate cell volumes
    cell_vol=zeros(ny,nx);
    % Calculate cell volumes based on the grid dimensions
    for i = 1:ny
        for j = 1:nx
            cell_vol(i,j) = dxVec(j) * dyVec(i); % Volume of each cell
        end
    end

   
    %Visualization of the grid.
    figure(1);
    % Plot vertical lines (constant j)
    plot(X, Y, 'Color', 'b', 'LineWidth', 0.5);
    hold on;
    % Plot horizontal lines (constant i)
    plot(X', Y', 'Color', 'b', 'LineWidth', 0.5);
    axis equal;
    xlabel('X');
    ylabel('Y');
    title(['Mesh for a flat plate with inflation layers at the' ...
        '     boundary and leading edge']);
   
   
else
    fprintf("DEFINE CORRECT VALUES FOR DOMAIN AND PLATE SIZES!!!")

end

%pause()

tic

%% VARIABLES AND COEFFITIENTS

%Main cell center pressure and velocity fields 
u=zeros(ny,nx); %Velocity in X axis
v=zeros(ny,nx);% Velocity in Y axis 
p=ones(ny,nx);%Pressure

%face velocities
u_face=zeros(ny,nx+1);
v_face=zeros(ny+1,nx);

u_face(:,1)=u0;

%Guessed nodal pressure and velocity fields
u_star=zeros(ny,nx); %X component
v_star=zeros(ny,nx); %Y component
%p_star=zeros(ny,nx); %Pressure


%Correction pressure  field
p_prime = zeros(ny,nx); %Pressure 

%Momentum equation coeffitients
aW=zeros(ny,nx); %Weast
aN=zeros(ny,nx); %North
aE=zeros(ny,nx); %East
aS=zeros(ny,nx); %South
aP=zeros(ny,nx); %Diagonal
aPv=zeros(1,nx);%Diagonal Coeffitienst for simetric boundary condition
%at the free-stream zones

%Pressure equation coeffitients
ap_W=zeros(ny,nx); %Weast
ap_N=zeros(ny,nx); %North
ap_E=zeros(ny,nx); %East
ap_S=zeros(ny,nx); %South
ap_P=zeros(ny,nx); %Diagonal

%Source terms 
suX=zeros(ny,nx); %Source for U
suY=zeros(ny,nx); %Source for V
suP=zeros(ny,nx); %Source for P

%Residual
rsid_x=zeros(ny,nx); %X momentum
rsid_y=zeros(ny,nx); %Y momentum
rsid_p=zeros(ny,nx);% Pressure correction eq.
rsid_cont=zeros(ny,nx); % Continuity

%error from residual
err_x=1;
err_y=1;
err_p=1;

%main iterations counter 
iterations_cont=0;

%Error for interior iterations
epsilon_uv=0.9e-15;
epsilon_p=5e-9;

%Underelaxation factors
alpha_uv=0.7;
alpha_p=0.005;

max_iterations=5000;% Max outer iterations <---------------------
max_iterations_uv=60;% Max iterations for momentum eqs.
max_iterations_p=120;% Max iterations for pressure eq.
residualsMat=zeros(max_iterations,3);% Residual Matrix
%iterations=0; %Iterations counter
error_tgt=1.139e-18; %Target error 5.1e-18 for Re=2000,1.2e-18
max_residual=1e10; % unstable values
convergedFlg=false; %Flag for convergence 


%Diffusive fluxes
dE=zeros(ny,nx);
for i=1:ny
    for j=1:nx
        if j==nx
            dE(i,j)=2*mu*dyVec(i)/dxVec(j);
        else
            epsilon_m=Xctrs(i,j+1)-Xctrs(i,j);
            dE(i,j)=mu*dyVec(i)/epsilon_m;
        end
    end
end

dW=zeros(ny,nx);

for i=1:ny
    for j=1:nx
        if j==1
            dW(i,j)=2*mu*dyVec(i)/dxVec(j);
        else
            epsilon_m=Xctrs(i,j)-Xctrs(i,j-1);
            dW(i,j)=mu*dyVec(i)/epsilon_m;
        end
    end
end

dN=zeros(ny,nx);

for i=1:ny
    for j=1:nx
        if i==1      
            dN(i,j)=2*mu*dxVec(j)/dyVec(i);
        else
            epsilon_n=Yctrs(i-1,j)-Yctrs(i,j);
            dN(i,j)=mu*dxVec(j)/epsilon_n;
        end
    end
end

dS=zeros(ny,nx);

for i=1:ny
    for j=1:nx
        if i==ny      
            dS(i,j)=2*mu*dxVec(j)/dyVec(i);
        else
            epsilon_s=Yctrs(i,j)-Yctrs(i+1,j);
            dS(i,j)=mu*dxVec(j)/epsilon_s;
        end
    end
end



%% SOLVER 

%nx_upstr,nx_fp,nx_dwnstr
%(j<=nx_upstr) || j>(nx-nx_dwnstr)

% Main loop
while convergedFlg==false
    %Iterations counter
    iterations_cont=iterations_cont+1;
    
    %Reset continuity residual for the new iteration
    rsid_cont(:) = 0;

    %1.- MOMENTUM LINK COEFFITIENTS AND SOURCES
    [aW,aE,aN,aS,aP,aPv,suX,suY] = momentum_link_coeff(nx,ny,...
    nx_upstr,nx_dwnstr,rho,dxVec,dyVec,u_face,v_face,p,u0,p0,...
    aW,aE,aN,aS,dW,dE,dN,dS,aP,aPv,suX,suY);

    %2.- SOLVE X MOMENTUM

    [u,rsid_x,err_x] = solve_momentum(err_x,max_iterations_uv,...
    nx,ny,u,alpha_uv,aP,aW,aN,aE,aS,suX,u_star,rsid_x);


    %3.- SOLVE Y MOMENTUM

    aP_temp_v=aP;% temporal central coeffitients matrix for v

    aP_temp_v(ny,nxSymcond)=aPv(nxSymcond);%Replace coeffitients that should not be
    %the same for the Simetry boundary condition

    [v,rsid_y,err_y] = solve_momentum(err_y,max_iterations_uv,...
    nx,ny,v,alpha_uv,aP_temp_v,aW,aN,aE,aS,suY,v_star,rsid_y);


    %4.- FACE VELOCITY COMPUTATION USING RIE - CHOW INTERPOLATION
    [u_face,v_face] = face_vel_intRC(nx,ny,nx_upstr,nx_dwnstr,u_face,...
        v_face,u,v,aP,aPv,alpha_uv,p,dxVec,dyVec,p0);


    %Neumman boundary condition for outlet flow at east edge
    u_face(:,end)=u(:,end);


    %5.- PRESSURE CORRECTION LINK COEFFITIENTS AND MASS INBALANCE (SOURCE)

    [ap_W,ap_N,ap_E,ap_S,ap_P,suP] = pressure_link_coeff(nx,ny,...
        nx_upstr,nx_dwnstr,dxVec,dyVec,u_face,v_face,aP,...
        aPv,ap_W,ap_N,ap_E,ap_S,ap_P,suP);

    
    %6.- SOLVE PRESSURE CORRECTION

    % Reset p_prime for the new iteration
    p_prime=zeros(ny,nx);

    [p_prime,rsid_p,err_p] = solve_presscorr(err_p,max_iterations_p...
        ,nx,ny,p_prime,ap_P,ap_W,ap_N,ap_E,ap_S,suP,rsid_p...
        ,epsilon_p);

    
    %7.- CORRECT PRESSURE

    p_star=p + alpha_p*p_prime;

    
    %8.- CORRECT CELL CENTER VELOCITY

    [u_star,v_star] = cvel_correct(nx,ny,nx_upstr,nx_dwnstr,dxVec,...
        dyVec,aP,aPv,u,v,p_prime);

        
    %9.- CORRECT FACE VELOCITY*

    [u_face,v_face] = fvel_correct(nx,ny,nx_upstr,nx_dwnstr,u_face,...
        v_face,alpha_uv,dxVec,dyVec,aP,aPv,p_prime);

    %10.- P=P_STAR
    p=p_star;



    %11.- CHECK RESIDUALS

    %Calculate residual continuity

    %for i=1:ny
    %    for j=1:nx
    %        rsid_cont(i,j)=dyVec(i)*(-u_face(i,j) + u_face(i,j+1)) +  ... 
    %            dxVec(j)*(v_face(i,j)-v_face(i+1,j));
    %    end
    %end

    %Calculate representative residual error continuity (sum of all cell 
    %mass balances)
    %rms(rsid_cont(:))
    %fprintf('%.20f',residualsMat(iterations_cont,3))
    
    continuity_scalF=0.000000070102575628;%0.000000070114960161;

    err_cont=abs(rms(suP(:)) - continuity_scalF);

    residualsMat(iterations_cont,:)=[err_x,err_y,err_cont];
    disp(residualsMat(iterations_cont,:))
    
    if err_x < error_tgt && err_y < error_tgt && err_cont < error_tgt
       convergedFlg=true;
       fprintf("Converged  at iteration\n")
       disp(iterations_cont)
       save('flatplateFields.mat', 'u', 'v' , 'p' );
    elseif err_x>max_residual || err_y > max_residual || err_cont > max_residual
        fprintf("Unstable solution iterations stopped at iteration ")
        disp(iterations_cont)
        break
    elseif iterations_cont >= max_iterations
        fprintf("Max iterations reached \n")
        disp(iterations_cont)
        break
    elseif isnan(err_x) || isnan(err_y) || isnan(err_cont)
        fprintf("The system became undetermined  at iteration \n")
        disp(iterations_cont)
        break
    else
        convergedFlg=false;
    end

end

%end loop

toc

%% POST PROCESS
fprintf('%s%d\n','Reynolds number ', Re)
fprintf('Total cells: %d\n', ncells);

res_init_point=1;% to skip high residual values at first iterations

figure(2)
plot(res_init_point:iterations_cont-1,...
    residualsMat(res_init_point:iterations_cont-1,1),...
    res_init_point:iterations_cont-1,...
    residualsMat(res_init_point:iterations_cont-1,2),...
    res_init_point:iterations_cont-1,...
    residualsMat(res_init_point:iterations_cont-1,3))
legend("U velocity","V velocity","Continuity ")
title("Residuals")
xlabel("Iterations")
ylabel("Residual")
yscale log


if convergedFlg==true
    
    %COMPUTATION OF OTHER VARIABLES OF INTEREST

    %VORTICITY
    [dvdx, dvdy] = gradient(v, Xctrs(1,:), Yctrs(:,1));
    [dudx, dudy] = gradient(u, Xctrs(1,:), Yctrs(:,1));
    vorticity = dvdx - dudy;

    %APROXIMATE STREAM FUNCTION

    psi = cumtrapz(Yctrs(:,1), u, 1); 

    %VELOCITY MAGNITUDE
    velocity_magnitude = sqrt(u.^2 + v.^2); 


    figure(3)
    contourf(Xctrs,Yctrs,u, 20, 'LineColor', 'none')
    title("Velocity in x axis (m/s)")
    xlabel("Height (m)")
    ylabel("Lenght (m)")
    colormap jet
    colorbar
    axis equal
    
    figure(4)
    contourf(Xctrs,Yctrs,v, 20, 'LineColor', 'none')
    title("Velocity in y axis (m/s)")
    xlabel("Height (m)")
    ylabel("Lenght (m)")
    colormap jet
    colorbar
    axis equal
    
    
    figure(5)
    contourf(Xctrs, Yctrs, velocity_magnitude, 20, 'LineColor', 'none')
    colorbar
    %hold on 
    %quiver(Xctrs,Yctrs,u,v,'k')
    title("Velocity vector field with magnitude m/s")
    xlabel("Height (m)")
    ylabel("Lenght (m)")
    colormap jet
    axis equal
    %hold off
    
    figure(6)
    contourf(Xctrs,Yctrs,p, 20, 'LineColor', 'none')
    title("pressure N/m^2")
    xlabel("Height (m)")
    ylabel("Lenght (m)")
    colormap jet
    colorbar
    axis equal

    figure(7)
    streamslice(Xctrs, Yctrs, u, v); % For a sparse representation
    title("streamlines")
    xlabel("Height (m)")
    ylabel("Lenght (m)")
    axis equal

    
    figure(8)
    contourf(Xctrs, Yctrs, vorticity, 20, 'LineColor', 'none');
    colorbar
    title('Vorticity Field')
    xlabel('X')
    ylabel('Y')
    axis equal
    colormap jet;
    

    figure(9)
    contour(Xctrs, Yctrs, psi, 30, 'LineWidth', 1.5);
    title('Streamfunction Lines')
    xlabel('X')
    ylabel('Y')
    axis equal

    figure(10)
    u_prof=[0,flip(u(:,nx_upstr + 5)'),u0];
    ycord_prof=[0,flip(yctrs),domHgt];
    plot(u_prof,ycord_prof,'- o')
    title("U profile of Flat Plate")
    xlabel("u velocity")
    ylabel("Y position")

    figure(11)
    contourf(Xctrs,Yctrs,cell_vol, 20, 'LineColor', 'none')
    title("Cell volumes")
    xlabel("Height (m)")
    ylabel("Lenght (m)")
    colormap jet
    colorbar
    axis equal

end

figure(12)
contourf(Xctrs,Yctrs,suP, 20, 'LineColor', 'none')
title("Continuity Residual (Kg*m^2/s)")
xlabel("Height (m)")
ylabel("Lenght (m)")
colorbar 
colormap jet
axis equal























