function check_model()
scope=setup_paths(); %#ok<NASGU>
[C,~]=study_config('parameter_BASE');
fprintf('Shared model: %s, N%d, dNuc=%g nm, c_fr CFL=%g\n',which('solver'),C.Nr,C.vofDNuc*1e9,C.cfl);
end
