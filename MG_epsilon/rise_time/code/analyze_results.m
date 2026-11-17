function summary=analyze_results()
scope=setup_paths(); %#ok<NASGU>
analyze_study('epsilon');
root=fileparts(fileparts(mfilename('fullpath')));summary=readtable(fullfile(root,'results','metrics','recomputed_summary.csv'));
end
