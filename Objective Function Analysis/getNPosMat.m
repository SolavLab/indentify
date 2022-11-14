function [time_mat, N_pos_mat,N_disp_mat] = getNPosMat(test)
% N_pos_mat(:,1,:) = (n_nodes)x(n_times) matrix of x positions
% each depth index k is [X,Y,Z] positions at time Tk
    if isfield(test,'pos_out')
        time_mat = test.pos_out.time;
        % disp_out at t=1 is used to restore pos_out at t=0, hence we need
        % to find to find the index correspondence between nodes in
        % disp_out and pos_out. Use only nodes that appear in both data
        % sets - 
        [~,ia,ib] = intersect(test.pos_out.ind,test.disp_out.ind); % test.pos_out.ind(ia)=test.disp_out.ind(ib)
        N_pos_mat = permute(cat(3,test.pos_out.x.data(ia,:),test.pos_out.y.data(ia,:),test.pos_out.z.data(ia,:)),[1 3 2]);
        N_disp_mat = permute(cat(3,test.disp_out.ux.data(ib,:),test.disp_out.uy.data(ib,:),test.disp_out.uz.data(ib,:)),[1 3 2]);
        % restore pos_out at t=0
        N_pos_mat(:,:,1) = N_pos_mat(:,:,2)-N_disp_mat(:,:,2);
        % delete data of all nodes that or not on the top and front
        % surfaces
        not_top_surface_nodes = find(N_pos_mat(:,2,1)~=0); %
        not_front_surface = find(N_pos_mat(:,3,1)~=0); %
        rows_to_delete = union(not_top_surface_nodes,not_front_surface);
        N_disp_mat(rows_to_delete,:,:) = [];
        N_pos_mat(rows_to_delete,:,:) = [];
    else
        error('pos_out is not a field of test');
    end  
    
    

    
% remove rows of indices that are not on the top front edge surface
%     to_delete_node_id = union(test.pos_out.ind(N_pos_mat(:,2,1)>0),test.pos_out.ind(abs(N_pos_mat(:,3,1))>1e-8));% indices of all nodes with y>0 or z~=0 at t=0
%     [~,~,pos_out_rows] = intersect(to_delete_node_id,test.disp_out.ind,'stable');    
%     N_pos_mat(to_delete_node_id,:,:) = [];
%     N_disp_mat_new(test.pos_out.ind(pos_out_rows),:,:) = [];
    
end