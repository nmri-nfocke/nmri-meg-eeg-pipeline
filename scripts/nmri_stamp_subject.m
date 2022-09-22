function [ subject ] = nmri_stamp_subject( subject, stamp_item, params )
%[ subject ] = nmri_stamp_subject( subject, stamp_item, params )
%   This will add a "stamp" to subject struct to mark a processing step has
%   been sucsessfully completed. If the same stamp is there, a warning will
%   be issued but overwritten


if ~exist('subject','var')
 error('Need subject info')
end

if ~exist('stamp_item','var')
 error('Need stamp info')
end

% check the params struct
if (~exist('params','var'))
 if (~exist('analysis_params.m','file'))
  error('Need to find analysis paramter file (analysis_params.m) in the current path, or have it in the call ')
 else
  analysis_params
  if (~exist('params','var')) 
   error('Problems with loading the paramter file (analysis_params.m)')  
  end
 end
end


% now check the subject struct
if (~isstruct(subject))
 % may be a mat file 
 if (exist(subject,'file'))
  if (strcmp(subject(end-1:end),'.m'))
   eval(subject);
  elseif (strcmp(subject(end-1:end),'.mat'))
   load(subject);
  else
   error('Subject provided is not a valid struct, .m or .mat file')
  end
 else
  error('Subject provided is not a valid struct, .m or .mat file')
 end
end

% make stamps branch if needed
if ~isfield(subject,'stamps')
 subject.stamps=[];
end


% make the stamp dir, if needed
if ~isfield(subject,'stamp_dir')
  subject.stamp_dir=fullfile(subject.analysis_dir,subject.id,'logs',subject.exam_id);
end
if (~exist(subject.stamp_dir,'dir'))
 mkdir(subject.stamp_dir)
end

% now check the stamp
if isfield(subject.stamps,stamp_item)
 fprintf(['Stamp = ' stamp_item ' already present - Will overwrite\n'])
end

subject.stamps.(stamp_item).params=params;
subject.stamps.(stamp_item).user=getenv('USER');
subject.stamps.(stamp_item).host=getenv('HOSTNAME');
% save versions of Matlab, Fieldtrip and SPM (if available)
subject.stamps.(stamp_item).version=['R' version('-release')];
[ftver, ftpath] = ft_version();
subject.stamps.(stamp_item).ft_version=ftver;
if ~isdeployed && exist('spm','file')
 spm_version=spm('ver');
 if strcmp(spm_version,'SPM12')
  spm_version=spm('Version');
 end
 subject.stamps.(stamp_item).spm_version=spm_version;
end
subject.stamps.(stamp_item).date=date;

% also save params valid during stamping
if exist('params','var')
 save(fullfile(subject.stamp_dir,stamp_item),'subject','params');
else
 save(fullfile(subject.stamp_dir,stamp_item),'subject');
end

fprintf('\nPlacing stamp %s for subject %s, %s\n\n',stamp_item,subject.id,subject.exam_id)

end

