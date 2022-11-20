function [test] = RunLog2Test(run_log)
% Initialize test{} struct cell array with metadata from
% run_log.mat/run_log.txt
    test = cell(numel(run_log.test),1); % Preallocating space
    for i=1:numel(run_log.test)
        test{run_log.test{i}.test_ind} = run_log.test{i};
    end
    test = loadDataFiles(test,1);
end



