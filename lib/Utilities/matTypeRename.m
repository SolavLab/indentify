% convert between old and new material model name formats
function [new_name] = matTypeRename(old_name)
    switch old_name
        case 'Ogden_symmetric'
            new_name = 'OM';
        case 'Ogden_1st'
            new_name = 'OG';
        otherwise
            new_name = old_name;
    end
end
          
