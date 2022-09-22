function [ all_metrics, all_layouts, all_subjects ] = nmri_load_metrics_sensor( all_subjects, metric_txt, freq, params )
% [ all_metrics, all_suma, all_subjects ] = nmri_load_metrics_sensor( all_subjects, metric_txt, freq, params )
%   will load a specific connectivity netric (or power/normalized power)
%   for a specific frequency band for all subjects given. Will skip missing
%   subjects (giving a warning)
%
%  will also do some basic processing of the metric (replace non-valid vertics wiht NaN, e.g. medial wall)
%
%  all_subjects  = cell array of subject structs
%  metric_txt    = char of metric (e.g. coh_img, or power, power_noise)
%  freq          = char of frequency or number (from analyis_params)
%  params        = optional, will read analyis_params if not given

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

if ~exist('metric_txt','var') || ~ischar(metric_txt)
 error('You need to specify which metric to use, e.g. ''coh_img'' for imaginary part of coherence')
end
if ~exist('freq','var') 
 error('You need to specify which frequency band to use, e.g. ''Alpha'' or the number based on the analysis_params')
end

% check freq text or number
if isnumeric(freq)
 freq_txt=params.freqsNames{freq};
 freq_n=freq;
elseif ischar(freq)
 freq_txt=freq;
 freq_n=find(strcmp(params.freqsNames,freq));
else
 error('Could not determine frequency requested')
end

if isstruct(all_subjects)
 % in case a single subject is given
 all_subjects={all_subjects};
end


all_metrics={};
all_layouts={};

for n=1:length(all_subjects)
 % parse over all subjects
 subject=all_subjects{n};
 fprintf('Subject=%s, Exam_ID=%s, Metric=%s, Freq=%s\n',subject.id,subject.exam_id,metric_txt,freq_txt)
 [ this_params ] = nmri_get_modality_params( params, subject.dtype );

 % in sensor space, we treat power as any other metric
 metric_file=fullfile(subject.stats_dir,['sensor_' metric_txt '_' subject.id '_' subject.exam_id '_' freq_txt '.mat']);
 if exist(metric_file,'file')
  load(metric_file,'metric','layout')
 else
  warning(sprintf('Could not load metric=%s for subject=%s, exam_id=%s, file=%s',metric_txt,subject.id,subject.exam_id,metric_file))
  all_subjects{n}=[];
  continue
 end
 

 % remove identy axis, if 2D
 if size(metric,2)>1 && size(metric,1)>1
  metric(eye(size(metric))==1)=NaN;
 end
 
 % and append
 all_metrics{end+1,1}=metric;
 all_layouts{end+1,1}=layout;
 
end
% now clean up skipped subjects
all_subjects(all(cellfun(@isempty,all_subjects),2),:)=[];
end

