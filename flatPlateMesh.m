%FLAT PLATE MESH 
%CODE  TO GENERATE MESH OVER A PLATE WITH INFLATION LAYERS AT THE 
%BOUNDARY, LEADING AND TRAILING EDGES 
%THERE IS A FREE - STREAM ZONE AT THE LEADING EDGE AND OTHER IN THE 
%TRAILING EDGE 

close all

%% GEOMETRICAL DATA 

%Geometry
platLgt=0.6;                    %Flat plate length (m)
domHgt=1;                       %Domain Height (m)
domLgt=1.8;                     %Domain Length (m)
distPlat=0.6;                   %Distance from origin to leading edge

if (distPlat+ platLgt)<domLgt
    %VALID GEOMETRY
    
    %NUMBER OF CELLS IN EACH SECTION
    nx=93;% number of cells in x axis 
    ny=45;% number of cells in y axis
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

    dy_1=0.00167;% Height of the first cell (m);0.00021,0.0001;%
    grthF=1.3;% Geometrical growth  Factor
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





