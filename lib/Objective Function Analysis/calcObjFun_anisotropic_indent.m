function [obj_val] = calcObjFun_anisotropic_indent(test,objectiveStruct)

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
force_sim = test.indenter_RB_out.Fz.data; % 1x(Nt)

%% Nodal Weights
% Calculate the distance between each node and the center of the indenter
if isfield(objectiveStruct,'nodeList')
    nodeList = objectiveStruct.nodeList;
    X_Wn = getWeight(pos_data,2,nodeList); %1xNodes
    Y_Wn = getWeight(pos_data,3,nodeList); %1xNodes
    Z_Wn = getWeight(pos_data,1,nodeList); %1xNodes
else
    n = size(pos_data, 1); % Determine the number of nodes
    defaultNodeList = true(1, n); % Default nodeList as a logical array of ones
    nodeList = defaultNodeList;
    X_Wn = getWeight(pos_data,2); %1xNodes
    Y_Wn = getWeight(pos_data,3); %1xNodes
    Z_Wn = getWeight(pos_data,1); %1xNodes
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

% Z displacement value
% weighed according to overall distance from indenter center
disp_res_mag = disp_exp_z-disp_sim_z(nodeList,:); % Not interpolated (synchronised+same mesh)
disp_exp_mag = max(abs(disp_exp_z(:,end))); %normalize by the greatest value overall
squared_normalized_disp_res = (disp_res_mag./disp_exp_mag).^2;
% Replace NaN values with 0 (normalizing by disp_exp_mag(:,:,1)==0)
squared_normalized_disp_res(isnan(squared_normalized_disp_res))=0;
node_weighted_squared_normalized_disp_res_Z = Z_Wn*squared_normalized_disp_res;
obj_val.Fu_z = node_weighted_squared_normalized_disp_res_Z; % Result: 1x(Nt)

% Y displacement value
% weighed according to distance from Y symmetry line
disp_res_mag = disp_exp_y-disp_sim_y(nodeList,:); % Not interpolated (synchronised+same mesh)
disp_exp_mag = max(abs(disp_exp_y(:,end))); %normalize by the greatest value overall
squared_normalized_disp_res = (disp_res_mag./disp_exp_mag).^2;
% Replace NaN values with 0 (normalizing by disp_exp_mag(:,:,1)==0)
squared_normalized_disp_res(isnan(squared_normalized_disp_res))=0;
node_weighted_squared_normalized_disp_res_Y = Y_Wn*squared_normalized_disp_res;
obj_val.Fu_y = node_weighted_squared_normalized_disp_res_Y;

% X displacement value
% weighed according to distance from X symmetry line
disp_res_mag = disp_exp_x-disp_sim_x(nodeList,:); % Not interpolated (synchronised+same mesh)
disp_exp_mag = max(abs(disp_exp_x(:,end))); %normalize by the greatest value overall
squared_normalized_disp_res = (disp_res_mag./disp_exp_mag).^2;
% Replace NaN values with 0 (normalizing by disp_exp_mag(:,:,1)==0)
squared_normalized_disp_res(isnan(squared_normalized_disp_res))=0;
node_weighted_squared_normalized_disp_res_X = X_Wn*squared_normalized_disp_res;
obj_val.Fu_x = node_weighted_squared_normalized_disp_res_X;


end