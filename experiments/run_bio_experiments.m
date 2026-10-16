function run_bio_experiments()
% Launch the author's original Bio1--Bio10 experiment entrypoint from repo root.
% This is a long-running experiment requiring MATLAB, COBRA, Gurobi, and MToP.
repo_root = fileparts(fileparts(mfilename('fullpath')));
assert(exist(fullfile(repo_root,'src','Bio_test','build_pool_bio.m'),'file')==2);
original_dir = pwd;
cleanup = onCleanup(@()cd(original_dir));
cd(repo_root);
run(fullfile(repo_root,'src','Bio_test','build_pool_bio.m'));
end
