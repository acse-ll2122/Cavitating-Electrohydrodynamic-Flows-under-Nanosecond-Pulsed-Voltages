function write_json(file,data)
    fid=fopen([file '.tmp'],'w');
    fprintf(fid,'%s',jsonencode(data));
    fclose(fid);
    movefile([file '.tmp'],file,'f');
end
