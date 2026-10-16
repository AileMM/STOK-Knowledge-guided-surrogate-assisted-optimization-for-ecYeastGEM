function verify_STOK_inputs()
% Structural/input checks; does not solve an FBA problem.
root = fileparts(fileparts(mfilename('fullpath')));
assert(isfile(fullfile(root,'models','kcats_inf.mat')),'Missing coefficient mapping.');
assert(isfile(fullfile(root,'models','ecYeastGEM_batch.mat')),'Missing base ecGEM model.');
k = load(fullfile(root,'models','kcats_inf.mat'),'kcats_inf');
assert(numel(k.kcats_inf)==4261,'Expected 4261 variable positions.');
for i=1:10
    kb = fullfile(root,'knowledge_base','benchmarks_bio',sprintf('bio%d-M2-d4261-k10.mat',i));
    assert(isfile(kb),'Missing Bio knowledge archive: %d',i);
    biofn = fullfile(root,'src','bio_objectives',sprintf('bio%d.m',i));
    assert(isfile(biofn),'Missing bio objective function: %d',i);
end
M = dir(fullfile(root,'models','bio_source_models','sourceModel_*.mat'));
assert(numel(M)==10,'Expected 10 environment-specific source models.');
fprintf('Structural checks passed: 4261 variables, 10 Bio objective files, 10 model files, 10 knowledge files.\n');
if exist('optimizeCbModel','file')~=2
    warning('STOK:COBRAMissing','COBRA Toolbox optimizeCbModel not on MATLAB path.');
end
if exist('UniformPoint','file')~=2
    warning('STOK:MToPMissing','MToP UniformPoint not on MATLAB path.');
end
end
