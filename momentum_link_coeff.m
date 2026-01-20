function [aW,aE,aN,aS,aP,aPv,suX,suY] = momentum_link_coeff(nx,ny,...
    nx_upstr,nx_dwnstr,rho,dxVec,dyVec,u_face,v_face,p,u0,p0,...
    aW,aE,aN,aS,dW,dE,dN,dS,aP,aPv,suX,suY)

%Coeffitients for momentum equations 
%   Convection + Diffusion + Sources Coeffitiens for the momentum equations
%in ordE(i,j)r to have them in this format considE(i,j)ring trasported variable phi:
% Ap*phi_p = Aw*phi_w + An*phi_n + Ae*phi_e + As*phi_s + Su 

%         BOUNDARY CONDITIONS FOR FLAT PLATE FLOW
%           DIRICHLET FIXED VELOCITY (U=U0 ,V=0)
%                                                      O
%I                                                     U
%N                                                     T
%L                                                     L
%E                                                     E   
%T                                                     T
% ___ SYMETRIC ___|___   NO SLIP ___|___ SYMETRIC  ____

%nx_upstr,nx_fp,nx_dwnstr


%1.- MOMENTUM LINK COEFFITIENTS
    
    %1.1- Interior Cells
    for i=2:ny-1
        for j=2:nx-1
    
            %Mass fluxes at faces 
            fW=-rho*dyVec(i)*u_face(i,j);
            fE=rho*dyVec(i)*u_face(i,j+1); 
            fN=rho*dxVec(j)*v_face(i,j); 
            fS=-rho*dxVec(j)*v_face(i+1,j); 
            
            %coeffitients using upwind squeme
            aW(i,j)=dW(i,j) + max(0,-fW);
            aE(i,j)=dE(i,j) + max(0,-fE);
            aN(i,j)=dN(i,j) + max(0,-fN);
            aS(i,j)=dS(i,j) + max(0,-fS);
            aP(i,j)=aW(i,j) + aE(i,j) + aN(i,j) + aS(i,j) +fW +fN +fE +fS;
    
            %Sources (Pressure dE(i,j)rivative)
            suX(i,j)=0.5*dyVec(i)*(p(i,j-1) -p(i,j+1));
            suY(i,j)=0.5*dxVec(j)*(p(i+1,j) -p(i-1,j));
    
        end
    end
    
    %1.2- Walls
    
    %1.2.1- Left wall - West (INLET)
    
    j=1;
    
    for i =2:ny-1
    
        %Mass fluxes at faces 
        fW=-rho*dyVec(i)*u_face(i,j); 
        fE=rho*dyVec(i)*u_face(i,j+1); 
        fN=rho*dxVec(j)*v_face(i,j); 
        fS=-rho*dxVec(j)*v_face(i+1,j);   
            
        %coeffitients using upwind squeme
        aW(i,j)=dW(i,j); 
        aE(i,j)=dE(i,j) + max(0,-fE);
        aN(i,j)=dN(i,j) + max(0,-fN);
        aS(i,j)=dS(i,j) + max(0,-fS);
        aP(i,j)=aW(i,j) + aE(i,j) + aN(i,j) + aS(i,j) +fN +fE +fS;
    
        %Sources (Pressure derivative)
        suX(i,j)=dyVec(i)*(p(i,j) - p(i,j+1)) + u0*(aW(i,j) -fW);
        suY(i,j)=0.5*dxVec(j)*(p(i+1,j) -p(i-1,j));
    
    end
    
    %1.2.2- Top wall - North  (DIRICHLET FIXED VELOCITY)
    
    i=1;
    
    for j=2:nx-1
    
        %Mass fluxes
        fW=-rho*dyVec(i)*u_face(i,j); 
        fE=rho*dyVec(i)*u_face(i,j+1); 
        %fN=rho*dxVec(j)*v_face(i,j); % No Mass fluxes at North face
        fS=-rho*dxVec(j)*v_face(i+1,j);   
            
        %coeffitients using upwind scheme
        aW(i,j)=dW(i,j) + max(0,-fW);
        aE(i,j)=dE(i,j) + max(0,-fE);
        aN(i,j)=dN(i,j);
        aS(i,j)=dS(i,j) + max(0,-fS);
        aP(i,j)=aW(i,j) + aE(i,j) + aN(i,j) + aS(i,j) + fW + fE + fS;
    
        %Sources (Pressure dE(i,j)rivative + boundary condition)
        suX(i,j)=0.5*dyVec(i)*(p(i,j-1) - p(i,j+1)) + u0*aN(i,j);
        suY(i,j)=0.5*dxVec(j)*(p(i+1,j) - p(i,j));
    
    end
    
    %1.2.3- Right Wall - East  (Outlet)
    
    j=nx;
    for i=2:ny-1
       
        %Mass fluxes at faces
        fW=-rho*dyVec(i)*u_face(i,j);  
        fE=rho*dyVec(i)*u_face(i,j+1); %Neumman condition 
        fN=rho*dxVec(j)*v_face(i,j); 
        fS=-rho*dxVec(j)*v_face(i+1,j);   
        
        %coeffitients using upwind squeme
        aW(i,j)=dW(i,j) + max(0,-fW);
        %aE(i,j)=0;
        aN(i,j)=dN(i,j) + max(0,-fN);
        aS(i,j)=dS(i,j) + max(0,-fS);
        aP(i,j)=aW(i,j) + aN(i,j) + aS(i,j)  +fW +fN +fS +fE;
        
        %Sources (Pressure dE(i,j)rivative)
        suX(i,j)=dyVec(i)*(0.5*(p(i,j) + p(i,j-1)) - p0);
        suY(i,j)=0.5*dxVec(j)*(p(i+1,j) -p(i-1,j));
    
    end
    
    %1.2.4- Bottom wall -- South (Symetric in the free-stream and no slip 
    % on the plate)
    
    i=ny;
    
    for j=2:nx-1

        if (j<=nx_upstr) || j>(nx-nx_dwnstr)

            %central coeffitient for momentum equation for v will be
            %different due to Symetric BC

            %Mass fluxes at faces 
            fW=-rho*dyVec(i)*u_face(i,j);
            fE=rho*dyVec(i)*u_face(i,j+1); 
            fN=rho*dxVec(j)*v_face(i,j); 
            %fS=-rho*dxVec(j)*v_face(i+1,j); No mass flux at south face

            %coeffitients using upwind squeme
            aW(i,j)=dW(i,j) + max(0,-fW);
            aE(i,j)=dE(i,j) + max(0,-fE);
            aN(i,j)=dN(i,j) + max(0,-fN);
            aS(i,j)=dS(i,j);%No appllicable for U
            aP(i,j)=aW(i,j) + aE(i,j) + aN(i,j) +fW +fN +fE;
            aPv(j)=aW(i,j) + aE(i,j) + aN(i,j)+ aS(i,j) +fW +fN +fE;

            %Sources (Pressure dE(i,j)rivative)
            suX(i,j)=0.5*dyVec(i)*(p(i,j-1) -p(i,j+1));
            suY(i,j)=0.5*dxVec(j)*(p(i,j) -p(i-1,j));

        else

            %central Coeffitient for u and v are the same
    
            %Mass fluxes at faces 
            fW=-rho*dyVec(i)*u_face(i,j);
            fE=rho*dyVec(i)*u_face(i,j+1); 
            fN=rho*dxVec(j)*v_face(i,j); 
            %fS=-rho*dxVec(j)*v_face(i+1,j); No mass flux at south face
            
            %coeffitients using upwind squeme
            aW(i,j)=dW(i,j) + max(0,-fW);
            aE(i,j)=dE(i,j) + max(0,-fE);
            aN(i,j)=dN(i,j) + max(0,-fN);
            aS(i,j)=dS(i,j);
            aP(i,j)=aW(i,j) + aE(i,j) + aN(i,j) + aS(i,j) +fW +fN +fE;
        
            %Sources (Pressure dE(i,j)rivative)
            suX(i,j)=0.5*dyVec(i)*(p(i,j-1) -p(i,j+1));
            suY(i,j)=0.5*dxVec(j)*(p(i,j) -p(i-1,j));

        end
    
    end
    
    %1.3- Corners
    
    %1.3.1- North-West
    i=1;
    j=1;
    
    %Mass fluxes at faces 
    fW=-rho*dyVec(i)*u_face(i,j); 
    fE=rho*dyVec(i)*u_face(i,j+1); 
    %fN=rho*dxVec(j)*v_face(i,j);  No mass flux at north face
    fS=-rho*dxVec(j)*v_face(i+1,j);
    
    %coeffitients using upwind squeme
    aW(i,j)=dW(i,j);
    aE(i,j)=dE(i,j) + max(0,-fE);
    aN(i,j)=dN(i,j);
    aS(i,j)=dS(i,j) + max(0,-fS);
    aP(i,j)=aW(i,j) + aE(i,j) + aN(i,j) + aS(i,j) + fE +fS;
    
    %Sources (Pressure dE(i,j)rivative + boundary condition)
    suX(i,j)=dyVec(i)*(p(i,j) - p(i,j+1)) + u0*(aW(i,j) + aN(i,j) -fW);
    suY(i,j)=0.5*dxVec(j)*(p(i+1,j) - p(i,j));
    
    %1.3.2- North-East
    i=1;
    j=nx;
    
    %Mass fluxes at faces 
    fW=-rho*dyVec(i)*u_face(i,j);
    fE=rho*dyVec(i)*u_face(i,j+1);
    %fN=rho*dxVec(j)*v_face(i,j);  No mass flux at north face
    fS=-rho*dxVec(j)*v_face(i+1,j); 
    
    %coeffitients using upwind squeme
    aW(i,j)=dW(i,j) + max(0,-fW);
    %aE(i,j)=2*dE(i,j);
    aN(i,j)=dN(i,j);
    aS(i,j)=dS(i,j) + max(0,-fS);
    aP(i,j)=aW(i,j) + aN(i,j) + aS(i,j) +fW +fS +fE;
    
    %Sources (Pressure dE(i,j)rivative + boundary condition)
    suX(i,j)=dyVec(i)*(0.5*(p(i,j) + p(i,j-1)) - p0) + u0*aN(i,j);
    suY(i,j)=0.5*dxVec(j)*(p(i+1,j) - p(i,j));
    
    %1.3.3- South-East (Symetric at South and outlet at East)
    i=ny;
    j=nx;
    
    %Mass fluxes at faces 
    fW=-rho*dyVec(i)*u_face(i,j);  
    fE=rho*dyVec(i)*u_face(i,j+1); 
    fN=rho*dxVec(j)*v_face(i,j); 
    %fS=-rho*dxVec(j)*v_face(i+1,j); No mass flux at south face
    
    %coeffitients using upwind scheme
    aW(i,j)=dW(i,j) + max(0,-fW);
    %aE(i,j)=2*dE(i,j);
    aN(i,j)=dN(i,j) + max(0,-fN);
    aS(i,j)=dS(i,j);%No applicable for u due to SBC
    aP(i,j)=aW(i,j) + aN(i,j) + fW + fN + fE;
    aPv(j)= aW(i,j) + aE(i,j) + aN(i,j) + aS(i,j) +fW +fN +fE;
    
    %Sources (Pressure dE(i,j)rivative)
    suX(i,j)=dyVec(i)*(0.5*(p(i,j) + p(i,j-1)) - p0);
    suY(i,j)=0.5*dxVec(j)*(p(i,j) -p(i-1,j));
    
    %1.3.4- South-West (Symetric at South and inlet at West)
    i=ny;
    j=1;
    
    %Mass fluxes at faces 
    fW=-rho*dyVec(i)*u_face(i,j);
    fE=rho*dyVec(i)*u_face(i,j+1); 
    fN=rho*dxVec(j)*v_face(i,j); 
    %fS=-rho*dxVec(j)*v_face(i+1,j);
    
    %coeffitients using upwind squeme
    aW(i,j)=dW(i,j);
    aE(i,j)=dE(i,j) + max(0,-fE);
    aN(i,j)=dN(i,j) + max(0,-fN);
    aS(i,j)=dS(i,j);%No applicable for u due to SBC
    aP(i,j)=aW(i,j) + aE(i,j) + aN(i,j)  +fN +fE;
    aPv(j)=aW(i,j) + aE(i,j) + aN(i,j) + aS(i,j)  +fN +fE;

    %Sources (Pressure dE(i,j)rivative)
    suX(i,j)=dyVec(i)*(p(i,j) - p(i,j+1)) + u0*(aW(i,j) - fW);
    suY(i,j)=0.5*dxVec(j)*(p(i,j) -p(i-1,j));
        
end

