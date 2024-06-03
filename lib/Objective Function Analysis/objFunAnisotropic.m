function [obj_val] = objFunAnisotropic(test,objectiveStruct)

%% get intenter's center of mass
if isfield(test,'MeshGeometry')
    indenter_center_of_mass_sim = test.MeshGeometry.Indenter.center_of_mass;
end
%% baseline ("experimental") data
pos_exp = objectiveStruct.pos_exp; % (Nnodes)x3x(Nt)
pos_exp_z = pos_exp.z.data(:,end);
pos_exp_y = pos_exp.y.data(:,end);
pos_exp_x = pos_exp.x.data(:,end);
force_exp = objectiveStruct.force_exp(:,end); % 1x(Nt)
%% trial ("simulated") data
pos_sim_z = test.pos_out.z.data(:,end);
pos_sim_y = test.pos_out.y.data(:,end);
pos_sim_x = test.pos_out.x.data(:,end);
force_sim = test.indenter_RB_out.Fz.data(:,end);
%% Nodal Weights
% Calculate the distance between each node and the center of the indenter
n = length(pos_exp_z);
distances = vecnorm([pos_exp_x';pos_exp_y';pos_exp_z']+repmat(indenter_center_of_mass_sim',1,n));
% Get the weight function corresponding to the distance - the further from
% the indenter the less weight this node will have on the result
[Wn] =getWn(distances,objectiveStruct.indenterRadius);
Wn = Wn/norm(Wn,2);
%% Residuals in force
force_res = force_exp-force_sim;
force_res_mag = abs(force_res);
force_exp_mag = abs(force_exp); % Normalization factor for force_res_mag
squared_normalized_force_res = (force_res_mag./force_exp_mag).^2;
% Replace NaN values with 0 (originating from normalizing by 0).
squared_normalized_force_res(isnan(squared_normalized_force_res))=0; % Ff(delta)@(p_sim,p_exp) eq. (3.1) in paper
%% Residuals in displacement Z
disp_res_mag = pos_exp_z-pos_sim_z; % Not interpolated (synchronised+same mesh)
disp_exp_mag = pos_exp_z;
squared_normalized_disp_res = (disp_res_mag./disp_exp_mag).^2;
% Replace NaN values with 0 (normalizing by disp_exp_mag(:,:,1)==0)
squared_normalized_disp_res(isnan(squared_normalized_disp_res))=0;
node_weighted_squared_normalized_disp_res_Z = (Wn)*squared_normalized_disp_res;
%% Residuals in displacement X
disp_res_mag = pos_exp_x-pos_sim_x; % Not interpolated (synchronised+same mesh)
disp_exp_mag = pos_exp_x;
squared_normalized_disp_res = (disp_res_mag./disp_exp_mag).^2;
% Replace NaN values with 0 (normalizing by disp_exp_mag(:,:,1)==0)
squared_normalized_disp_res(isnan(squared_normalized_disp_res))=0;
node_weighted_squared_normalized_disp_res_X = (Wn)*squared_normalized_disp_res;
node_weighted_squared_normalized_disp_res_X2 = (1./(pos_exp_y+0.1))'*squared_normalized_disp_res;
%% Residuals in displacement Y
disp_res_mag = pos_exp_y-pos_sim_y; % Not interpolated (synchronised+same mesh)
disp_exp_mag = pos_exp_y;
squared_normalized_disp_res = (disp_res_mag./disp_exp_mag).^2;
% Replace NaN values with 0 (normalizing by disp_exp_mag(:,:,1)==0)
squared_normalized_disp_res(isnan(squared_normalized_disp_res))=0;
node_weighted_squared_normalized_disp_res_Y = (Wn)*squared_normalized_disp_res;
node_weighted_squared_normalized_disp_res_Y2 = (1./(pos_exp_x+0.1))'*squared_normalized_disp_res;
%% Set output structure
obj_val.Ff = squared_normalized_force_res(end);
obj_val.Fu_z = node_weighted_squared_normalized_disp_res_Z;
obj_val.Fu_x = node_weighted_squared_normalized_disp_res_X;
obj_val.Fu_y = node_weighted_squared_normalized_disp_res_Y;
obj_val.Fu_x2 = node_weighted_squared_normalized_disp_res_X2;
obj_val.Fu_y2 = node_weighted_squared_normalized_disp_res_Y2;

end
