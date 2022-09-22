function [ subject ] = nmri_read_subject(subj_id,exam_id,analysis_dir,raw_dir,raw_dataset)
%subject=nmri_read_subject(subj_id,exam_id,analysis_dir,raw_dir,raw_dataset)
%  
% This function will read a MEG or EEG dataset from the
% raw_dir/subj_id/exam_id directory and create the necessary directories,
% links and.m files (if not present) in analysis_dir. This is the basis
% for further processing with nmri_... functions
%
% subj_id   =   Subject ID, usually basename
% exam_id   =   realtive path to exam, e.g. meg/<subj_id>.MEG.rest.ds
%               will check for path without <subj_id> as well
% analysis_dir = root dir for analyis, is not set will take pwd
% raw_dir   =   where to search for raw files, will try /data/storage and 
%               pwd, if not set

% written by NF 10/2016 - 07/2018

% some presets, may be editable later

% check the call
if (nargin<2) 
    error('Need at least subj_id and exam_id')
end
if ( ~exist('analysis_dir','var') || isempty(analysis_dir) ) 
    analysis_dir=pwd;
end

if ( ~exist('raw_dir','var') ) 
   raw_dir=pwd;
end


% check if we have a subject_info m.m file
if (exist(fullfile(analysis_dir,subj_id,'subject_info.m'),'file'))
 disp('Found subject_info file - loading');
 clear subject
 %save path
 curpath=pwd;
 cd(fullfile(analysis_dir,subj_id))
 eval('subject_info');
 cd(curpath);
else
 % make empty struct
 subject=[];
end

if ( exist('raw_dataset','var') ) 
   subject.raw_dataset=raw_dataset;
end

% deal with paths
paths=[];
if ~isempty(getenv('NMRI_STORAGE'))
 paths.storage=getenv('NMRI_STORAGE');
else
 paths.storage='/data/storage';
end




% read params file
[ params ] = nmri_get_params( analysis_dir );

% check for the raw file

% check if set in our current subject struct
if (~isfield(subject,'raw_dataset') || ~exist(subject.raw_dataset,'file'))
 % search for raw dataset
 if (exist(fullfile(raw_dir,subj_id,exam_id),'file'))
  subject.raw_dataset=fullfile(raw_dir,subj_id,exam_id);
 elseif (exist(fullfile(analysis_dir,subj_id,exam_id),'file'))
  subject.raw_dataset=fullfile(analysis_dir,subj_id,exam_id);  
 elseif (exist(fullfile(paths.storage,subj_id,exam_id),'file'))
  subject.raw_dataset=fullfile(paths.storage,subj_id,exam_id); 
 else    
  % now check with basename again
  [epath, exam]=strtok(exam_id,'/');
  if (~isempty(exam))
   exam=exam(2:end);
  end
  if (exist(fullfile(raw_dir,subj_id,epath,[subj_id '.' exam]),'file'))
   subject.raw_dataset=fullfile(raw_dir,subj_id,epath,[subj_id '.' exam]);
   exam_id=fullfile(epath,[subj_id '.' exam]);
  elseif (exist(fullfile(analysis_dir,subj_id,epath,[subj_id '.' exam]),'file'))
   subject.raw_dataset=fullfile(analysis_dir,subj_id,epath,[subj_id '.' exam]);   
   exam_id=fullfile(epath,[subj_id '.' exam]);
  elseif (exist(fullfile(paths.storage,subj_id,epath,[subj_id '.' exam]),'file'))
   subject.raw_dataset=fullfile(paths.storage,subj_id,epath,[subj_id '.' exam]);  
   exam_id=fullfile(epath,[subj_id '.' exam]);  
  elseif (exist(fullfile(raw_dir),'file'))
   subject.raw_dataset=fullfile(raw_dir);
  else
   error(['No raw dataset found for this subject'])
  end
 end
end


% deal with gziped edf
if strcmpi(subject.raw_dataset(end-6:end),'.edf.gz')

 [status,flags]=fileattrib(subject.raw_dataset);
 if status && flags.UserWrite
  % unzip in place if write permission
  fprintf('Unzippig in place=%s\n',subject.raw_dataset);
  system(['gzip -d ''' subject.raw_dataset ''''])
  if exist(subject.raw_dataset(1:end-3),'file')
   subject.raw_dataset=subject.raw_dataset(1:end-3);
  else
   error('Unzip of EDF failed, maybe file corrupt / incomplete?')
  end
 else
  % cannot unzip in place, so copy to analyis dir
  [pa fi ext]=fileparts(subject.raw_dataset);
  new_dir=fullfile(analysis_dir,subj_id,'eeg');
  if ~exist(new_dir,'dir')
    mkdirp(new_dir);
  end
  new_target=fullfile(new_dir,fi);
  fprintf('Unzippig in analysis dir=%s\n',new_target);
  system(['gzip -dc ' subject.raw_dataset ' > ' new_target])
  if exist(new_target,'file')
   subject.raw_dataset=new_target;
  else
   error('Unzip of EDF failed, maybe file corrupt / incomplete?')
  end
 end
end


% now deal with the special nuiscance of EGI's JAR reader...they always
% want to open the files with random access, i.e. need write permissions

if strcmpi(subject.raw_dataset(end-3:end),'.mff')
 mff_reader=nmri_check_mff_reader(subject.raw_dataset);
 fprintf('Using MFF-Reader=%s\n',mff_reader);
 if ~strcmp(mff_reader,'egi_mff_v1')
  % the VAR Reader needs write permissions, check if we have that
  [status,flags]=fileattrib(subject.raw_dataset);
  if ~status || ~flags.UserWrite
   % copy to analyis dir
   new_target=fullfile(analysis_dir,subj_id,exam_id);
   if ~exist(new_target,'dir')
    mkdirp(new_target);
    system(['cp -r ' subject.raw_dataset '/* ' new_target])
   end
   fprintf('Copying MFF to analysis dir=%s\n',new_target)
   if exist(new_target,'dir')
    [status,flags]=fileattrib(new_target);
    if status && flags.UserWrite
     system(['chmod -R u+w ' new_target]);
    end
    subject.raw_dataset=new_target;
   else
    error('Could not copy .mff file to analysis dir - probably permission problem')
   end
  end
 end
end



% put additional info in subject
subject.id=subj_id;
%clean exam_id
exam_id=strrep(exam_id,'.ds','');
exam_id=strrep(exam_id,'.mff','');
exam_id=strrep(exam_id,'.fif','');
exam_id=strrep(exam_id,'.eeg','');
exam_id=strrep(exam_id,'.edf','');
exam_id=strrep(exam_id,'.EDF','');
exam_id=strrep(exam_id,'.gz','');
exam_id=strrep(exam_id,'/','_');
exam_id=strrep(exam_id,'\','_');
exam_id=strrep(exam_id,'(','_');
exam_id=strrep(exam_id,')','');
exam_id=strrep(exam_id,[subj_id '.'],'');
exam_id=strrep(exam_id,[subj_id '_'],'');
exam_id=strrep(exam_id,[subj_id],'');
exam_id=strrep(exam_id,'.','_');
exam_id=strrep(exam_id,'__','_');
subject.exam_id=exam_id;
subject.analysis_dir=analysis_dir;



% call central dataset selector
subject=nmri_determine_datatype(subject);

% we need dtype as a minimum
if (~isfield(subject,'dtype') || isempty(subject.dtype))
 error('Could not detect the datatype of this data. Please check the fileformat / file access. Will need to stop here.')
end


% read hdr, if not there
if (~isfield(subject,'hdr') || isempty(subject.hdr))
 error('HDR data not found, should have been read by nmri_determine_datatype before, investigate')
end


if (~isfield(subject,'valid_channels'))
 switch subject.dtype
 case 'MEG'
  subject.valid_channels={'MEG'};
 case 'EEG'
  subject.valid_channels={'EEG'};
 case 'EEG_invasive'
  subject.valid_channels={'EEG'};
 otherwise
  error('No supported datatype detected')
 end 
end

% make conf dir
if (~exist(fullfile(analysis_dir,subj_id,'conf',subject.exam_id),'dir'))
 mkdir(fullfile(analysis_dir,subj_id,'conf',subject.exam_id))
end

% save subject as JSON , if not present
if (~exist(fullfile(analysis_dir,subj_id,'conf',subject.exam_id,'subject_info.json'),'file'))
 subject_json=subject;
 if isfield(subject_json,'layout') && isstruct(subject_json.layout)
  % not intended for JSON edit, too long
  subject_json=rmfield(subject_json,'layout');
 end
 if isfield(subject_json,'hdr') && isstruct(subject_json.hdr)
  % not intended for JSON edit 
  subject_json=rmfield(subject_json,'hdr');
 end
 nf_json_write(fullfile(subject.analysis_dir,subject.id,'conf',subject.exam_id,'subject_info.json'),subject_json);
end

% save subject file, if not present
if (~exist(fullfile(analysis_dir,subj_id,'conf',subject.exam_id,'subject_info.m'),'file'))

 
 fid=fopen(fullfile(analysis_dir,subj_id,'conf',subject.exam_id,'subject_info.m'),'w');
 fprintf(fid,'%% This is an automatically generated subject file\n');
 fprintf(fid,'%% you should modify it if needed (bad channels, trial_funtion,...)\n\n');
 fprintf(fid,'\n\n%%NOTE: This file is for posteriority only.\n%%Pipelines after 2019 DO NOT READ THIS any more\n\n');
 
 fprintf(fid,'subject=[];\n');
 fprintf(fid,'subject.raw_dataset=''%s'';\n',subject.raw_dataset);
 fprintf(fid,'subject.id=''%s'';\n',subject.id);
 fprintf(fid,'subject.analysis_dir=''%s'';\n',subject.analysis_dir);
 fprintf(fid,'subject.exam_id=''%s'';\n',subject.exam_id);
 
 switch subject.dtype
 case 'MEG'
  fprintf(fid,'subject.dtype=''%s'';\n',subject.dtype);
  fprintf(fid,'subject.valid_channels={''MEG''}; %% use all MEG channels\n');
  subject.valid_channels={'MEG'};
  fprintf(fid,'%%subject.valid_channels={''MEG'',''-MRO22''}; %% exclude certain channels (in this case MRO22)\n');
 case 'EEG'
  fprintf(fid,'subject.dtype=''%s'';\n',subject.dtype);
  fprintf(fid,'subject.valid_channels={''EEG''}; %% use all EEG channels\n');
  subject.valid_channels={'EEG'};
  fprintf(fid,'%%subject.valid_channels={''EEG'',''-Cz''}; %% exclude certain channels (in this case Cz)\n');
 case 'EEG_invasive'
  fprintf(fid,'subject.dtype=''%s'';\n',subject.dtype);
  fprintf(fid,'subject.valid_channels={''EEG''}; %% use all EEG channels\n');
  subject.valid_channels={'EEG'};
  fprintf(fid,'%%subject.valid_channels={''EEG'',''-LTM1''}; %% exclude certain channels (in this case LTM1)\n');
 end 

 if isfield(subject,'layout') && ischar(subject.layout)
  fprintf(fid,'subject.layout=''%s''; %% Channel layout file\n',subject.layout);
 end
 if isfield(subject,'montages') && ischar(subject.layout)
  fprintf(fid,'subject.montages=''%s''; %% Montages file for viewer\n',subject.montages);
 end
 
 fprintf(fid,'%%subject.trial_fun=''''; %% if a specific trial function should be used\n');
 fprintf(fid,'%%subject.bad_trials={}; %% list trials that are known to be bad\n');

 fprintf(fid,'\n%% Further comments for this Subject can be added below\n');
 fprintf(fid,'%%\n');
 fclose(fid);
end

% also save a .mat file

save(fullfile(analysis_dir,subj_id,'conf',subject.exam_id,'subject_info.mat'),'subject')





    

end

