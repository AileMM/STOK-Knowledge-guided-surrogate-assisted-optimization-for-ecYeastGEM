function startup_STOK(mtop_root)
% Prepare MATLAB search path for STOK and its external MToP dependency.
% Usage: startup_STOK('C:\path\to\MTO-Platform')
repo_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(repo_root,'src')));
addpath(fullfile(repo_root,'models'));
addpath(fullfile(repo_root,'models','bio_source_models'));
addpath(fullfile(repo_root,'knowledge_base','benchmarks_bio'));
addpath(fullfile(repo_root,'knowledge_base','benchmarks_generic'));
addpath(fullfile(repo_root,'experiments'));
if nargin>=1 && ~isempty(mtop_root)
    assert(isfolder(mtop_root),'STOK:MToPPath','MToP directory does not exist');
    if isfolder(fullfile(mtop_root,'MTO'))
        addpath(genpath(fullfile(mtop_root,'MTO')),'-end');
    else
        addpath(genpath(mtop_root),'-end');
    end
else
    warning('STOK:MToPMissing','MToP v1.5 path not supplied. Utilities such as UniformPoint and W1_Optimum may be missing.');
end
fprintf('STOK MATLAB paths prepared: %s\n',repo_root);
end
