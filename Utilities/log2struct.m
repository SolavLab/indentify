function [s] = log2struct(varargin)
    switch nargin
        case 1 
           [time, data, data_label] =  importFEBio_logfile(varargin{1});
        case 2
           [time, data, data_label] =  importFEBio_logfile(varargin{1},varargin{2});
        case 3
           [time, data, data_label] =  importFEBio_logfile(varargin{1},varargin{2},varargin{3});
    end
    data_label_arr =strsplit(data_label,';');
    N_labels = numel(data_label_arr);
    s = struct;
    s.N = size(data,1);
    s.time = time;
    if ~(varargin{3}) % adds index column to data
        s.ind = data(:,1,end);
        data(:,1,:) = []; % remove index column from data
    end
    for i=1:N_labels      
        s.(data_label_arr{i}).data = squeeze(data(:,i,:)); % row_ind->elements, columns->time
        if strcmp(data_label_arr{i},'J')
                s.(data_label_arr{i}).data(:,1) = 1; % initialize J=1 @t=0
        end
        if size( s.(data_label_arr{i}).data,2)==1
             s.(data_label_arr{i}).data =  s.(data_label_arr{i}).data';
        end
    end
    
    
end

