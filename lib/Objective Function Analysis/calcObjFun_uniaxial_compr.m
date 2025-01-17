function [obj_val] = calcObjFun_uniaxial_compr(test,objectiveStruct, useForceOnly)

if useForceOnly
    force_exp = objectiveStruct.force_exp; % 1x(Nt)
    force_sim = abs(sum(test.force_out.Rz.data,1)); % 1x(Nt)
    % Force value
    force_res = force_exp-force_sim;
    force_res_mag = abs(force_res);
    force_exp_mag = max(abs(force_exp)); % Normalization factor for force_res_mag
    squared_normalized_force_res = (force_res_mag./force_exp_mag).^2;
    % Replace NaN values with 0 (originating from normalizing by 0).
    squared_normalized_force_res(isnan(squared_normalized_force_res))=0;
    obj_val.Ff = squared_normalized_force_res;
    obj_val.Fu_r = 0;

else
    pos_data = objectiveStruct.pos_data; % (Nnodes)x3

    %% baseline ("experimental") data

    disp_exp = objectiveStruct.disp_exp; % (Nnodes)x3
    disp_exp_z = disp_exp.uz.data;
    disp_exp_y = disp_exp.uy.data;
    disp_exp_x = disp_exp.ux.data;
    force_exp = objectiveStruct.force_exp; % 1x(Nt)

    %% trial ("simulated") data
    % pos_sim = test.pos_out; % (Nnodes)x3

    disp_sim = test.disp_out;       % (Nnodes)x3x(Nt)
    disp_sim_z = disp_sim.uz.data;
    disp_sim_y = disp_sim.uy.data;
    disp_sim_x = disp_sim.ux.data;
    force_sim = sum(test.force_out.Rz.data,1); % 1x(Nt)

    %% Nodal Weights
    % Calculate the distance between each node and the center of the indenter
    if isfield(objectiveStruct,'nodeList')
        nodeList = objectiveStruct.nodeList;
    else
        n = size(pos_data, 1); % Determine the number of nodes
        defaultNodeList = true(1, n); % Default nodeList as a logical array of ones
        nodeList = defaultNodeList;
    end

    %% Values
    % This section takes all the data and returns a 1x(Nt) vector of the
    % objective function values at each time step.

    % Force value
    force_res = force_exp-force_sim;
    force_res_mag = abs(force_res);
    force_exp_mag = max(abs(force_exp)); % Normalization factor for force_res_mag
    squared_normalized_force_res = (force_res_mag./force_exp_mag).^2;
    % Replace NaN values with 0 (originating from normalizing by 0).
    squared_normalized_force_res(isnan(squared_normalized_force_res))=0;
    obj_val.Ff = squared_normalized_force_res;

    % Radial displacement value
    % weighed according to overall distance from indenter center
    radial_disp_exp = sqrt(disp_exp_x.^2+disp_exp_y.^2);
    radial_disp_sim = sqrt(disp_sim_x.^2+disp_sim_y.^2);
    disp_res_mag = radial_disp_exp-radial_disp_sim(nodeList,:); % Not interpolated (synchronised+same mesh)
    disp_exp_mag = max(abs(radial_disp_exp(:,end))); %normalize by the greatest value overall
    squared_normalized_disp_res = (disp_res_mag./disp_exp_mag).^2;
    % Replace NaN values with 0 (normalizing by disp_exp_mag(:,:,1)==0)
    squared_normalized_disp_res(isnan(squared_normalized_disp_res))=0;
    obj_val.Fu_r = mean(squared_normalized_disp_res,1); % Result: 1x(Nt)
end

end