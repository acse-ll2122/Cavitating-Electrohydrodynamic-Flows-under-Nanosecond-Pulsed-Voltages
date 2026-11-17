function G=parameter_groups(root)
% Layout membership is separate from catalog physical values.
G=jsondecode(fileread(fullfile(root,'results','display_groups.json')));
T=readtable(catalog_file(),'TextType','string');T=T(T.kind=="parameter",:);
columns={'Cc','Cv','d_nm','n0_m3','pth_MPa'};
for j=1:numel(G)
 for k=1:numel(G(j).ids),G(j).levels(k)=T.(columns{j})(T.id==string(G(j).ids{k}));end
end
end
