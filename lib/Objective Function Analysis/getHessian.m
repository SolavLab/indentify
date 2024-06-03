function [H] = getHessian(varargin)

% function [H]=getHessian(X,Z,n,target_ind,leap)
% ------------------------------------------------------------------------

% X is an n dimensional cell array where each element is a vector
% indicating the parameter values in space

% Z is an n dimensional array of same size as X where each element is the
% objective function value for the set of parameters indexed by X.

% n: truncation error order O(h^n) (up to 8)

%target_ind: this is the position of the parameter point within the objective function space

% leap: distances between adjacent data points for calculation differences

%% Parse input 

switch nargin
    case 2
        X=varargin{1};
        Z=varargin{2};
        target_ind=ceil(size(X)/2); %pick center point in space
        n=2; %default truncation error
        leap=1; %default leap
    case 3
        X=varargin{1};
        Z=varargin{2};
        n=varargin{3};
        target_ind=ceil(size(X)/2);
        leap=1;
    case 4
        X=varargin{1};
        Z=varargin{2};
        n=varargin{3};
        target_ind=varargin{4};
        leap=1;
    case 5
        X=varargin{1};
        Z=varargin{2};
        n=varargin{3};
        target_ind=varargin{4};
        leap=varargin{5};
end

%% Compute Hessian

[~,~,~,~,W_xy]=getXYCoeffMat(n); %Mixed numeric stencil coefficients (eq. (7) in Thesis)
Wtt = secondDerWeights(n); %Second derivative stencil coefficients
r = n/2; %the amount of indexes from either side of the target point
H = zeros(length(target_ind)); %initialize Hessian matrix
v = 1:ndims(X); %define vector of dimensions in the numeric space

for ii = 1:length(target_ind)
    for jj = ii:length(target_ind)
        try
            if ii == jj %diagonal values of the hessian require the second derivative in that direction
                C=num2cell(target_ind);
                C{ii}=':';
                local_Z = Z(C{:}); %get 1D slice!
                local_X = X(C{:}); %get 1D slice!
                p_i = target_ind(ii);
                target = cell2mat(local_X(p_i));
                step = cell2mat(local_X(p_i+leap));
                h_pi = norm(step-target);
                drv_values=local_Z((p_i-r*leap):leap:(p_i+r*leap));
                dev=0;
                for q=1:length(drv_values)
                    dev=dev+drv_values(q)*Wtt(q);
                end
                H(ii,ii) = dev/(h_pi^2);
            else
                operateOn = [ii,jj]; %choose the dimensions that the local hessian will operate on
                temp_target = target_ind; temp_target(operateOn)=[];
                temp_target = [target_ind(operateOn),temp_target]; %set the indices such that the first two correlate to d^2/(dp_i)(dp_j)
                local_X = permute(X, [operateOn, setdiff(v, operateOn)]); %set the array such that the first and second dimensions define a matrix to operate on
                local_Z = permute(Z, [operateOn, setdiff(v, operateOn)]);
                C=num2cell(temp_target);
                C{1}=':'; C{2}=':'; %structure the input local_X(:,:,target_ind(other dimensions))
                local_X = local_X(C{:}); % get 2D slice!
                local_Z = local_Z(C{:}); % get 2D slice!
                % get numeric step
                p_i = temp_target(1);
                p_j = temp_target(2);
                target = cell2mat(local_X(p_i,p_j)); %get vector of values from x
                step_pj = cell2mat(local_X(p_i,p_j+leap));
                step_pi = cell2mat(local_X(p_i+leap,p_j));
                h_pj = norm(step_pj-target); % delta_pj
                h_pi = norm(step_pi-target); % delta_pi
                %calculate the hessian!
                H(ii,jj) = sum((W_xy.*local_Z((p_i-r*leap):leap:(p_i+r*leap),(p_j-r*leap):leap:(p_j+r*leap))), 'all')/(h_pi*h_pj);
                H(jj,ii) = H(ii,jj);
            end
        catch ME % optimal parameter set is on the boundary of the parameters range
            H = [];
            warning('Optimal parameter set is on the boundary of the parameters range');
        end
    end
end