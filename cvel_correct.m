function [u_star,v_star] = cvel_correct(nx,ny,nx_upstr,nx_dwnstr,...
    dxVec,dyVec,aP,aPv,u,v,p_prime)
%CORRECTION OF CELL CENTER VELOCITIES 
% using pressure corection field to correct pressure velocities 

u_star=zeros(ny,nx);
v_star=zeros(ny,nx);

%Velocity in X axis 

    %Interior cells 
    for i=1:ny
        for j=2:nx-1
            u_star(i,j)=u(i,j) + 0.5*(dyVec(i)/aP(i,j))*(p_prime(i,j-1) - ... 
                p_prime(i,j+1));
        end
    end
    
    %Left
    j=1;
    for i=1:ny
        u_star(i,j)=u(i,j) + 0.5*(dyVec(i)/aP(i,j))*(p_prime(i,j)-p_prime(i,j+1));
    end
    
    %Right
    j=nx;
    for i=1:ny
        u_star(i,j)=u(i,j) + 0.5*(dyVec(i)/aP(i,j))*(p_prime(i,j-1)-p_prime(i,j));
    end
    
    %Velocity in Y axis
    
    %interior cells
    for i=2:ny-1
        for j=1:nx
            v_star(i,j)= v(i,j) + 0.5*(dxVec(j)/aP(i,j))*(p_prime(i+1,j) - ...
                p_prime(i-1,j)); 
        end
    end
    
    %Top
    
    i=1;
    for j=1:nx
        v_star(i,j)=v(i,j) + 0.5*(dxVec(j)/aP(i,j))*(p_prime(i+1,j) - p_prime(i,j));
    end
    
    %Bottom
    i=ny;
    for j=1:nx
        
        if (j<=nx_upstr) || j>(nx-nx_dwnstr)% Symetric BC  coeffitients
            v_star(i,j)=v(i,j) + 0.5*(dxVec(j)/aPv(j))*(p_prime(i,j) - p_prime(i-1,j));
        else       
            v_star(i,j)=v(i,j) + 0.5*(dxVec(j)/aP(j))*(p_prime(i,j) - p_prime(i-1,j));
        end
    end

end