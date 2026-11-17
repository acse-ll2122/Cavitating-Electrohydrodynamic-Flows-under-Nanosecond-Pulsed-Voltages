function T=comparison_cases(kind)
% The 50 kV / 1 ns pair is shared by the two comparison groups.
T=readtable(catalog_file(),'TextType','string');T=T(T.kind=="epsilon",:);
if strcmp(kind,'voltage'),T=T(T.rise_ns==1,:);
elseif strcmp(kind,'rise_time'),T=T(T.U_kV==50,:);
else,error('EHD:Comparison','Use voltage or rise_time.');end
end
