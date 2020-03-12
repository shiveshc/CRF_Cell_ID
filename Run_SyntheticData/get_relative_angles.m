%%% function to calculate angle between pair of neurons in data and
%%% pair of neurons in atlas. This angle defines edge potential.
function angle_matrix = get_relative_angles(X_rot,Y_rot,Z_rot,X,Y,Z,node1,node2)
    
    p_prime = [X(node1,1),Y(node1,1),Z(node1,1)];
    q_prime = [X(node2,1),Y(node2,1),Z(node2,1)];
    
    angle_matrix = zeros(size(X_rot,1),size(X_rot,1));
    for m = 1:size(X_rot,1)
        for n = m:size(X_rot,1)
            if n == m
                angle_matrix(m,n) = -1;
            else
                p = [X_rot(m,1),Y_rot(m,1),Z_rot(m,1)];
                q = [X_rot(n,1),Y_rot(n,1),Z_rot(n,1)];
                angle_matrix(m,n) = (p-q)*(p_prime - q_prime)'/(norm(p-q)*norm(p_prime-q_prime));
                angle_matrix(n,m) = -angle_matrix(m,n);
            end
        end
    end
%     angle_matrix = (1+angle_matrix)/2;
end
        