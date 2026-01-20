function [u_face,v_face] = fvel_correct(nx,ny,nx_upstr,nx_dwnstr,u_face,...
    v_face,alpha_uv,dxVec,dyVec,aP,aPv,p_prime)
%FACE VELOCITY CORRECTION
%face velocity correction

    %Velocity in X axis 
    for i =1:ny
        for j=2:nx
            u_face(i,j)=u_face(i,j) +...
                0.5*alpha_uv*dyVec(i)*(1/aP(i,j-1) + ...
                1/aP(i,j))*(p_prime(i,j-1)-p_prime(i,j));
        end
    end
    
    %Velocity in Y axis 
    for i=2:ny
        for j=1:nx


            if ((j<=nx_upstr) || j>(nx-nx_dwnstr)) && (i==ny)
                v_face(i,j)=v_face(i,j) + 0.5*alpha_uv*dxVec(j)*(1/aP(i-1,j) + ...
                1/aPv(j))*(p_prime(i,j)-p_prime(i-1,j));
            else
                v_face(i,j)=v_face(i,j) + 0.5*alpha_uv*dxVec(j)*(1/aP(i-1,j) + ...
                1/aP(i,j))*(p_prime(i,j)-p_prime(i-1,j));
            end

            
        end
    end

end