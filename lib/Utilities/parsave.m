function parsave(fname, varname)
    temp_struct.(varname) = evalin('caller',varname);
    save(fname,'-struct','temp_struct')
end