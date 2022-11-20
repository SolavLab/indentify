function [test] = getLogfileNames(test,febio_spec)
    if isfield(febio_spec,'Output')  % Add output file names to test structure
        if isfield(febio_spec.Output,'logfile')
            if isfield(febio_spec.Output.logfile, 'node_data')
                for ii=1:numel(febio_spec.Output.logfile.node_data)
                    temp_name = febio_spec.Output.logfile.node_data{ii}.ATTR.file;
                    test.node_data_files{ii} = temp_name;
                end
            end
            if isfield(febio_spec.Output.logfile, 'element_data')
                for ii=1:numel(febio_spec.Output.logfile.element_data)
                    temp_name = febio_spec.Output.logfile.element_data{ii}.ATTR.file;
                    test.element_data_files{ii} = temp_name;
                end
            end
            if isfield(febio_spec.Output.logfile, 'rigid_body_data')
                for ii=1:numel(febio_spec.Output.logfile.rigid_body_data)
                    temp_name = febio_spec.Output.logfile.rigid_body_data{ii}.ATTR.file;
                    test.rigid_body_data_files{ii} = temp_name;
                end
            end
        end
    end
end