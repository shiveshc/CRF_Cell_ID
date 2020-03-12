%%% function to calculate angle between pair of neurons in data and
%%% pair of neurons in atlas. This angle defines edge potential. This
%%% function is based on updated angle relationships based on annotated
%%% data

function angle_matrix = get_relative_angles_updated_atlas(X_rot,Y_rot,Z_rot,X,Y,Z,node1,node2,angle_vec_atlas)
    
    p_prime = [X(node1,1),Y(node1,1),Z(node1,1)];
    q_prime = [X(node2,1),Y(node2,1),Z(node2,1)];
    
    angle_matrix = angle_vec_atlas*((p_prime - q_prime)'/norm(p_prime-q_prime));
    angle_matrix = reshape(angle_matrix,size(X_rot,1),size(X_rot,1));
    
    for m = 1:size(angle_matrix,1)
        angle_matrix(m,m) = -1;
    end
end
        
        