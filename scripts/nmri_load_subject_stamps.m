function [ subject ] = nmri_load_subject_stamps( subject, params )
%[ subject ] = nmri_load_subject_stamps( subject , params )
%   Will load subject struct from the stamps/log directory and unified
%   cache in subjec conf

% check time
runT=now;

% check if subject struct is okay
if ~isfield(subject,'id') || ~isfield(subject,'exam_id') || ~isfield(subject,'analysis_dir') || ~isfield(subject,'dtype')
 subject_status = false;
 warning(['Subject struct is lacking basic info -- should not happen'])
else
 % this seems to be at least basically correct
 
 % check the params struct
 if (~exist('params','var'))
  if isfield(subject,'params')
   params=subject.params;
  else
   prev_dir=pwd;
   if exist(subject.analysis_dir,'dir')
    cd(subject.analysis_dir)
   else
    warning(['Could not CD into the original analyis dir (' subject.analysis_dir ') to load the analysis params. Stay where we are...'])
   end
   if (~exist('analysis_params.m','file'))
    error('Need to find analysis paramter file (analysis_params.m) in analysis dir ')
   else
    analysis_params
    if (~exist('params','var')) 
     error('Problems with loading the paramter file (analysis_params.m)')  
    end
   end
   cd(prev_dir)
  end
 end
 [ params ] = nmri_get_modality_params( params, subject.dtype );

 % check if the current analyis dir make sense, i.e. still exists
 if ~exist(fullfile(subject.analysis_dir,subject.id,'conf',subject.exam_id,'subject_info.mat'),'file')
  % if not, try in pwd instead
  if ~exist(fullfile(pwd,subject.id,'conf',subject.exam_id,'subject_info.mat'),'file')
   error('Could not access the subject_info.mat file in the specified dir and the pwd. Cannot work with ')
  else
   subject_file=fullfile(pwd,subject.id,'conf',subject.exam_id,'subject_info.mat');
  end
 else
  subject_file=fullfile(subject.analysis_dir,subject.id,'conf',subject.exam_id,'subject_info.mat');
 end
 
 % preserve pre-load analyis dir
 analysis_dir=subject.analysis_dir;


  
 % make sure we have the stamp dir
 if ~isfield(subject,'stamp_dir')
  subject.stamp_dir=fullfile(subject.analysis_dir,subject.id,'logs',subject.exam_id);
 end
 if (~isfield(subject,'stamps'))
  subject.stamps=[];
 end
 
 
 % check the subjects master cache

 subject_cache=fullfile(subject.analysis_dir,subject.id,'conf',subject.exam_id,'subject_cache.mat');
 
 
 if exist(subject_cache,'file')
  cached_subject=load(subject_file);
  % for now, empty timestamp
  cached_subject.max_timestamp=0;  
  if ~isfield(cached_subject,'subject')
   % should not happen, unless file is broken somehow
   % make an empty struct
   cached_subject.subject=[];
   cached_subject.timestamp.basic=0;
   cached_subject.max_timestamp=0;  
   cached_subject.source=[];
  end
  if ~isfield(cached_subject,'timestamp')
   cached_subject.timestamp.basic=0; 
  end
  if ~isfield(cached_subject,'max_timestamp')
   cached_subject.max_timestamp=0; 
  end
  if ~isfield(cached_subject,'source')
   cached_subject.source=[];
  end
 else
  % make an empty struct
  cached_subject.subject=subject;
  cached_subject.timestamp.basic=0;
  cached_subject.max_timestamp=0;  
  cached_subject.source=[];
 end

 
 % now check all the available stamps
 %fprintf('Reading stamps dir...\n')
 %tic
 allS=dir(fullfile(subject.stamp_dir,'*.mat'));
 %toc
 % sort by datenum, ascending is the default
 %[~,idx]=sort([allS.datenum]);
 
 % extract the real stamp name without extension
 for i=1:length(allS)
  [~,fi,~]=fileparts(allS(i).name);
  allS(i).stamp=fi;
 end
 
 % now check first in a "classical" way 
 items={'preproc','readingevents','artifactrejection','artifactrejectionICA','MRI_ctf','MRI_seg_SUMA','vigilance','hdmSUMA_nmriSUMA150_fs6_openmeeg','processing_sensor'};
 if isfield(params,'headmodel')
  items=[items {['hdmSUMA_' params.headmodel] ['hdmSUMA_individual_' params.headmodel],['processing_' params.headmodel],['processing_individual_' params.headmodel]}];
 end
 
 % check the intersection
 if ~isempty(allS)
  present=intersect(items,{allS.stamp},'stable'); % ordered as in items 
  % and invert, so last is read first
  present=flip(present);
 else
  % in case no stamps a present, we cannot use the stamp field
  present={};
 end
 stamps=[];
 
 % the classical items will be parsed by the logic of the analyis flow, not
 % the date
 for i=1:length(present)
  idx=find(strcmp(present{i},{allS.stamp}));
  % check if this file is newer than our current max timestamp
  if allS(idx).datenum > cached_subject.max_timestamp
  % load the file, this seems newer
   %fprintf('Loading stamps %s...\n',allS(idx).name)
   %tic
   load_file=fullfile(allS(idx).folder,allS(idx).name);
   tmp=load(load_file,'subject');
   % re-root
   if ~strcmp(tmp.subject.analysis_dir,analysis_dir)
    tmp.subject=nmri_reroot_subject(tmp.subject,analysis_dir);
   end
   % remember loaded from
   tmp.subject.loaded_from=load_file;
   %toc
   % get rid of data sometimes stored in the hdr for MFF
   if isfield(tmp.subject,'hdr') && isfield(tmp.subject.hdr,'orig') && isfield(tmp.subject.hdr.orig,'data')
    tmp.subject.hdr.orig=rmfield(tmp.subject.hdr.orig,'data');
   end
   if isfield(tmp.subject,'hdr') && isfield(tmp.subject.hdr,'orig') && isfield(tmp.subject.hdr.orig,'event')
    tmp.subject.hdr.orig=rmfield(tmp.subject.hdr.orig,'event');
   end
   % deal with stamps seperately, avoid recursion
   if isfield(tmp.subject,'stamps')
    if isempty(stamps)
     % this is the first stamp, so remember
     stamps=tmp.subject.stamps;
    end
    % and then remove
    tmp.subject=rmfield(tmp.subject,'stamps');
   end
   % remember for stamps
   if ~isfield(stamps,allS(idx).stamp)
    stamps.(allS(idx).stamp)=tmp.subject;
   end

   % and update
   %tic
   %fprintf('Updating cache...\n')
   fields=fieldnames(tmp.subject);
   for ii=1:length(fields)
    if ~isfield(cached_subject.subject,fields{ii}) || ~isequal(cached_subject.subject.(fields{ii}),tmp.subject.(fields{ii}))
     % not present or different, so updated
     cached_subject.subject.(fields{ii})=tmp.subject.(fields{ii});
     cached_subject.source.(fields{ii})=allS(idx).stamp;
     cached_subject.timestamp.(fields{ii})=allS(idx).datenum;
    end
   end
   % now make a new max timestamp
   fields=fieldnames(cached_subject.timestamp);
   tp=zeros(1,length(fields));
   for i=1:length(fields)
    tp(i)=getfield(cached_subject.timestamp,fields{i});
   end
   cached_subject.max_timestamp=max(tp);
   %toc 
  end
 end
 
 % now add the additional / canonical / non-standard stamps
 if isfield(allS,'stamp')
  other_stamps=setdiff({allS.stamp},fieldnames(stamps));
  for i=1:length(other_stamps)
   idx=find(strcmp(other_stamps{i},{allS.stamp}));
   %fprintf('Loading and adding stamp %s...\n',allS(idx).name)
   %tic
   tmp=load(fullfile(allS(idx).folder,allS(idx).name),'subject');
   % get rid of data sometimes stored in the hdr for MFF
   if isfield(tmp.subject,'hdr') && isfield(tmp.subject.hdr,'orig') && isfield(tmp.subject.hdr.orig,'data')
    tmp.subject.hdr.orig=rmfield(tmp.subject.hdr.orig,'data');
   end
   if isfield(tmp.subject,'hdr') && isfield(tmp.subject.hdr,'orig') && isfield(tmp.subject.hdr.orig,'event')
    tmp.subject.hdr.orig=rmfield(tmp.subject.hdr.orig,'event');
   end
   % avoid recursion
   if isfield(tmp.subject,'stamps') 
    tmp.subject=rmfield(tmp.subject,'stamps');
   end
   % add this stamp
   stamps.(allS(idx).stamp)=tmp.subject;
   %toc
  end
 
  % and add
  cached_subject.subject.stamps=stamps;
 end

 % deal with reroot
 if (~strcmp(cached_subject.subject.analysis_dir,analysis_dir))
  cached_subject.subject=nmri_reroot_subject(cached_subject.subject,analysis_dir); 
 end

 % get rid of data sometimes stored in the hdr for MFF
 if isfield(cached_subject.subject,'hdr') && isfield(cached_subject.subject.hdr,'orig') && isfield(cached_subject.subject.hdr.orig,'data')
  cached_subject.subject.hdr.orig=rmfield(cached_subject.subject.hdr.orig,'data');
 end
 if isfield(cached_subject.subject,'hdr') && isfield(cached_subject.subject.hdr,'orig') && isfield(cached_subject.subject.hdr.orig,'event')
  cached_subject.subject.hdr.orig=rmfield(cached_subject.subject.hdr.orig,'event');
 end
  
 
 % now make sure we have all hdm classed loaded (for new style hdm)
 if ~isfield(cached_subject.subject,'hdm_classes')
  % make sure we have the field
  cached_subject.subject.hdm_classes=[];
 end
 if isfield(cached_subject.subject,'stamps') && ~isempty(cached_subject.subject.stamps)
  f=fieldnames(cached_subject.subject.stamps);
  for ii=1:length(f)
   if isfield(cached_subject.subject.stamps.(f{ii}),'hdm_classes') && ~isempty(cached_subject.subject.stamps.(f{ii}).hdm_classes)
    hf=fieldnames(cached_subject.subject.stamps.(f{ii}).hdm_classes);
    for iii=1:length(hf)
     % now we add or update the hdm_class from the stamp
     if ~isfield(cached_subject.subject.hdm_classes,hf{iii})
      cached_subject.subject.hdm_classes.(hf{iii})=cached_subject.subject.stamps.(f{ii}).hdm_classes.(hf{iii});
     else
      % update and merge with existing
      cached_subject.subject.hdm_classes.(hf{iii})=nf_update_struct(cached_subject.subject.hdm_classes.(hf{iii}),cached_subject.subject.stamps.(f{ii}).hdm_classes.(hf{iii}));
     end
    end
   end
  end
 end
 
 
 % now make a unified timestamp
 fields=fieldnames(cached_subject.timestamp);
 tp=zeros(1,length(fields));
 for i=1:length(fields)
  tp(i)=getfield(cached_subject.timestamp,fields{i});
 end
 cached_subject.max_timestamp=max(tp);
 
 
 
 % now save in subject_cache.mat file
 try
  save(subject_cache,'-struct','cached_subject');
 catch
  fprintf('Could not write subject cache file in %s, probably a permission problem',subject_cache)
 end
 

 % now make active
 subject=cached_subject.subject;
 
 
 
 % make sure we always have some relevant fields
 if (~isfield(subject,'mri_dir'))
  subject.mri_dir=fullfile(subject.analysis_dir,subject.id,'mri');
 end
 if (~isfield(subject,'mri_ctf'))
  subject.mri_ctf=fullfile(subject.mri_dir,['mri_ctf_' subject.id '.mat']);
 end
 if (~isfield(subject,'mri_seg'))
  subject.mri_seg=fullfile(subject.mri_dir,['mri_seg_' subject.id '.mat']);
 end
 if (~isfield(subject,'suma_surface'))
  subject.suma_surface=fullfile(subject.mri_dir,['suma_surface_' subject.id '.mat']);
 end
 if (~isfield(subject,'clean_dataset'))
  subject.clean_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['clean_' subject.id '_' subject.exam_id '.mat']);
 end
 if (isfield(params,'useICA_clean') && params.useICA_clean==1)
  if (~isfield(subject,'cleanICA_dataset'))
   subject.cleanICA_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['cleanICA_' subject.id '_' subject.exam_id '.mat']);
  end
 end
 if (~isfield(subject,'stats_dir'))
  subject.stats_dir=fullfile(subject.analysis_dir,subject.id,'stats');
 end
 if (~isfield(subject,'stats') && isfield(params,'headmodel'))
  subject.stats=fullfile(subject.stats_dir,['source_stats_' subject.id '_' subject.exam_id '_' params.headmodel '.mat']);
 end
 if (~isfield(subject,'sensor_stats') && isfield(params,'headmodel'))
  subject.sensor_stats=fullfile(subject.stats_dir,['sensor_stats_' subject.id '_' subject.exam_id '_' params.headmodel '.mat']);
 end
 
end

fprintf('Total runtime for %s was %s (min:sec.ms)\n',subject.id,datestr(now-runT,'MM:SS.FFF'))
end

