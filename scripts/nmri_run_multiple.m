function [ SGEdirs]= nmri_run_multiple(subjects,cfg,cmd,params)
% [ SGEdirs , sge_ids]= nmri_run_multiple(subjects,cfg,cmd,params)
%   Function to run jobs via nmri_parrun & nmri_qsub on multiple subjects
%
%   Sets some defaults for MEG/EEG processing and makes sure params are
%   attached to subject structs
%
%   will do some basic checking, if job is still running (via lockfiles in
%   SGEdir).



if ~iscell(subjects)
 error('Subjects list needs to be a cell array')
end

if ~iscell(cmd)
 error('Command list needs to be a cell array, to be passed to nmri_qsub')
end


% make the config
if ~exist('cfg','var')
 cfg=[];
end
if ~isfield(cfg,'compile')
 % default is to compile
 cfg.compile=1;
end
if ~isfield(cfg,'fieldtrip')
 % default is to use fieldtrip
 cfg.fieldtrip=1;
end
if ~isfield(cfg,'overwrite')
 % default is not to overwrite
 cfg.overwrite=0;
end
if ~isfield(cfg,'rootdir')
 cfg.rootdir=pwd;
end

if ~isfield(cfg,'run_time')
 % default runtime is 2h
 cfg.run_time=120;
end

if ~isfield(cfg,'base_title')
 % default base title is cmd joined
 cfg.base_title=strjoin(cmd,'_');
end

if ~isfield(cfg,'chunck_size')
 % default to have one job / chunk
 cfg.chunck_size=1;
end


N=length(subjects);
SGEdirs=cell(N,1);

cfg.spm=1;
% now loop over all subjects
for n=1:N
 subject=subjects{n};
 % now make the individual title
 cfg.title=[cfg.base_title '_' subject.id '_' subject.exam_id];

 % check if running
 cfg.SGEdir=fullfile(cfg.rootdir,'SGE_calls',cfg.base_title,[subject.id '_' subject.exam_id]);

 if exist(fullfile(cfg.SGEdir,'run_lock'),'file') || exist(fullfile(cfg.SGEdir,'queue_lock'),'file')
  % found a lock file
  fprintf('Found a probably running job for ID=%s in %s, skipping this subject\n',subject.id,cfg.SGEdir)
 else
  % run it
  % attach params
  subject.params=params;
  % make paths absolute
  subject=nmri_path_absolute(subject);
  [SGEdirs{n},~,cfg_out]=nmri_qsub(cfg,cmd,subject);
  % remember the compiled hash to make sure this is constant and faster
  cfg.compiled_hash=cfg_out.compiled_hash;
 end

end

end
