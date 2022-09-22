function [ all_subjects ] = nmri_export_metrics_sensor( all_subjects, opt, params )
%[ all_subjects ] = nmri_export_metrics_sensor( all_subjects, opt, params )
% Will save any number of metrics / frequency bands (after scaling / 
% dim_reduction/ abs) in a 4D .mgh file. The order will be as in the 
% all_subject cell-array of structs
%
% this is the analogon to nmri_export_metrics, but runs in sensor space


% opt               = struct array of options
%  .layout          = fieldtrip layout structure to use as a reference, if
%                     not set, will take from 1st subject
%  .metrics         = cell array of metrics to use
%                     will take from params (+ power, power_noise), if not set
%  .freqs           = cell array of frequencies to use (either char name or index)
%                     will take all from params, if not set
%  .dim_reduction   = cell array of statisitcal measures to use for
%                     dimension reduction
%                     default, {'mean'}, i.e. mean 
%                     available: mean, median, none
%  .global          = cell array of global value modification
%                     'abs': absolute, 'pos': only postive, 'neg': only
%                     negative, '' or 'none':nothing
%  .scale           = cell array of global scaling
%                     'mean', 'median', 'rms': root-mean-square, 
%                     '' or 'none':nothing
%  .interpolate     = true / false, interpolate missing / NaN values 
%  .save_freesurfer = true / false, conform to freesurfer convention
%                     i.e. lh/rh and only cortex, default: false
%  .save_average    = true / false, write average (for all specified 
%                     save options), default: false
%  .save_mgh        = save as .mgh files - for Freesurfer / PALM
%                     default: true
%  .save_mat        = save as .mat files - for Matlab 
%                     default: true
%  .save_grp_mean   = save the mean per map/variation as .mgh
%                     default: false
%
%  .output          = output dir, will place 4D .mgh files there 
%                     default: <analysis_dir>/export/<date>
%  

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


% check opt and set defaults
if ~exist('opt','var')
 % defaults
 opt=[];
end




if ~isfield(opt,'metrics') || ~iscell(opt.metrics) || isempty(opt.metrics)
 % as a default we take everything
 opt.metrics={};
 all_metrics=[params.con_method {'power','power_noise'}];
 for i=1:length(all_metrics)
  if strcmp(all_metrics{i},'coh')
   opt.metrics{end+1}=[all_metrics{i} '_real'];
   opt.metrics{end+1}=[all_metrics{i} '_img'];
  else
   opt.metrics{end+1}=all_metrics{i};
  end
 end
end

if ~isfield(opt,'layout') 
 % load from 1st subject, assuming this to be the same for all
 opt.layout = [];
end



if ~isfield(opt,'save_freesurfer') || ~islogical(opt.save_freesurfer ) || isempty(opt.save_freesurfer )
 opt.save_freesurfer=false;
end
if ~isfield(opt,'save_average') || ~islogical(opt.save_average ) || isempty(opt.save_average )
 opt.save_average=false;
end
if ~isfield(opt,'save_mgh') || ~islogical(opt.save_mgh ) || isempty(opt.save_mgh )
 opt.save_mgh=true;
end
if ~isfield(opt,'save_mat') || ~islogical(opt.save_mat ) || isempty(opt.save_mat )
 opt.save_mat=true;
end
if ~isfield(opt,'save_grp_mean') || ~islogical(opt.save_grp_mean ) || isempty(opt.save_grp_mean )
 opt.save_grp_mean=false;
end

if ~isfield(opt,'dim_reduction') || ~iscell(opt.dim_reduction) || isempty(opt.dim_reduction)
 opt.dim_reduction={'mean'};
end

if ~isfield(opt,'scale') || ~iscell(opt.scale) || isempty(opt.scale)
 opt.scale={'none','mean'};
end

if ~isfield(opt,'global') || ~iscell(opt.global) || isempty(opt.global)
 opt.global={'abs','none'};
end

if ~isfield(opt,'interpolate') || ~islogical(opt.interpolate ) || isempty(opt.interpolate )
 opt.interpolate=false;
end

if ~isfield(opt,'output') || ~ischar(opt.output) || isempty(opt.output)
 opt.output='<analysis_dir>/export/<date>';
end

if ~isfield(opt,'freqs') || ~iscell(opt.freqs) || isempty(opt.freqs)
 % as a default we take everything
 opt.freqs=1:length(params.freqs);
end
% check freq text or number
freq_txt=cell(length(opt.freqs),1);
freq_n=zeros(length(opt.freqs),1);
for i=1:length(opt.freqs)
 if isnumeric(opt.freqs(i))
  freq_txt(i)=params.freqsNames(opt.freqs(i));
  freq_n(i)=opt.freqs(i);
 elseif iscell(opt.freqs(i)) && ischar(opt.freqs{i})
  freq_txt{i}=opt.freqs{i};
  freq_n=find(strcmp(params.freqsNames,opt.freqs{i}));
 else
  error('Could not determine frequency requested')
 end
end

if isstruct(all_subjects)
 % in case a single subject is given
 all_subjects={all_subjects};
end



% set for some basic info
subject=all_subjects{1};

% make outpur dir
outputdir=opt.output;
outputdir=strrep(outputdir,'<analysis_dir>',subject.analysis_dir);
outputdir=strrep(outputdir,'<date>',date);
if ~exist(outputdir,'dir')
 mkdir(outputdir)
end

% Main Loop - Per Metric and Per Frequency
for m=1:length(opt.metrics)
 for f=1:length(opt.freqs)
  fprintf('Metric=%s, Freq=%s\n',opt.metrics{m},freq_txt{f})

   % Load the data - for all subjects
  [all_metrics, all_layouts, subjects_back]=nmri_load_metrics_sensor(all_subjects,opt.metrics{m},freq_n(f),params);
  if ~isequal(subjects_back,all_subjects)
   warning(['Could not find the metric ' opt.metrics{m} '/' freq_txt{f} ' for all subjects -- skipping'])
   continue
  end
  
  % parse the layouts
  if isempty(opt.layout)
   % take the 1st as default
   opt.layout=all_layouts{1};
  end
  % load the layout
  if ischar(opt.layout)
   l=load(opt.layout);
   opt.layout=l.layout;
  end
  
  
  for si=1:length(all_layouts)
   % map this subject's metrics according to the master layout
   layout_mapping=nan(size(all_metrics{si},1),1);
   for li=1:length(all_layouts{si}.label)
    idx=find(strcmpi(all_layouts{si}.label{li},opt.layout.label));
    if length(idx)==1
     layout_mapping(li)=idx;
    else
     layout_mapping(li)=NaN;
     fprintf(['Could not find a unique mapping from channel (' all_layouts{si}.label{li} ') in the data and the layout.\n'])
    end
   end
   
   % then remap to layout mapping
   metric=all_metrics{si};
   % 2D
   if size(metric,2)>1 && size(metric,1)>1
    mapped_metric=nan(length(opt.layout.label));
    mapped_metric(layout_mapping(~isnan(layout_mapping)),layout_mapping(~isnan(layout_mapping)))=metric(~isnan(layout_mapping),~isnan(layout_mapping));
    all_metrics{si}=mapped_metric;
   elseif size(metric,2)==1 || size(metric,1)==1
    % 1D
    mapped_metric=nan(length(opt.layout.label),1);
    mapped_metric(layout_mapping(~isnan(layout_mapping)))=metric(~isnan(layout_mapping));
    all_metrics{si}=mapped_metric;
   else
    error('Could not identify dim size of metric')
   end
  end
  
  % start with loop over global abs/pos/neg
  for global_i=1:length(opt.global)
   % now apply global 
   fprintf('Global=%s\n',opt.global{global_i})
   global_metrics=all_metrics;
   for mi=1:length(global_metrics)
    switch(opt.global{global_i})
    case 'abs'
     global_metrics{mi}=abs(global_metrics{mi});
    case 'pos'
     global_metrics{mi}(global_metrics{mi}<0)=NaN;
    case 'neg'
     global_metrics{mi}(global_metrics{mi}>0)=NaN;
    end
   end
   
   % now loop over all dim_reductions
   for dim_i=1:length(opt.dim_reduction)
    % now check and dim reduce each metric
    fprintf('Dim Reduction=%s\n',opt.dim_reduction{dim_i})
    dim_metrics=global_metrics;
    for mi=1:length(dim_metrics)
     dim_metrics{mi}=reduce_dim(dim_metrics{mi},opt.dim_reduction{dim_i});
     % and rotate if needed
     if (size(dim_metrics{mi},2)>size(dim_metrics{mi},1))
      dim_metrics{mi}=dim_metrics{mi}';
     end
    end 

    % now loop over scale
    for scale_i=1:length(opt.scale)
     fprintf('Scale=%s\n',opt.scale{scale_i})
     scale_metrics=dim_metrics;
     for mi=1:length(scale_metrics)
      switch(opt.scale{scale_i})
       case 'mean'
        val=nanmean(reshape(scale_metrics{mi},numel(scale_metrics{mi}),1));
       case 'median'
        val=median(reshape(scale_metrics{mi},numel(scale_metrics{mi}),1),'omitnan');
       case 'rms'
        val=rms(reshape(scale_metrics{mi},numel(scale_metrics{mi}),1));
       otherwise
        val=1;
      end
      scale_metrics{mi}=scale_metrics{mi}/val;
     end
     
     % interpolate missing values if asked for
     if opt.interpolate
      icfg=[];
      icfg.layout=opt.layout;
      for mi=1:length(scale_metrics)
       [scale_metrics{mi},icfg]=nmri_interpolate_map(icfg,scale_metrics{mi});
      end
     end
    
     % prepare export
     if strcmp(opt.global{global_i},'pos') 
      global_txt='pos_only';
     elseif strcmp(opt.global{global_i},'neg')
      global_txt='neg_only';
     elseif strcmp(opt.global{global_i},'abs')
      global_txt='abs';
     else
      global_txt='no_abs';
     end
     if strcmp(opt.scale{scale_i},'mean') 
      scale_txt='mean_scaled';
     elseif strcmp(opt.scale{scale_i},'median')
      scale_txt='median_scaled';
     elseif strcmp(opt.scale{scale_i},'rms')
      scale_txt='rms_scaled';
     else
      scale_txt='not_scaled';
     end
     if strcmp(opt.dim_reduction{dim_i},'mean') 
      dim_txt='';
     elseif strcmp(opt.dim_reduction{dim_i},'none') 
      dim_txt='_fulldim'; 
     else
      dim_txt=['_dimred_' opt.dim_reduction{dim_i}];
     end
     
     if opt.interpolate
      interp_txt='_interp';
     else
      interp_txt='';
     end

     filebase=fullfile(outputdir,[ opt.metrics{m} '_' freq_txt{f} '_' global_txt '_' scale_txt dim_txt interp_txt '_N' num2str(length(all_subjects))]);
     if opt.save_mat
      fprintf('Saving as .mat...')
      save([filebase '.mat'],'scale_metrics','all_subjects')
      fprintf('done\n')
     end
     
     if opt.save_mgh && ~strcmp(opt.dim_reduction{dim_i},'none') 
      % now save as raw MGH  - as is
      fprintf('Saving as full 4D MGH/suma_all...')
      nmri_write_mgh([filebase '_all.mgh'],eye(4),scale_metrics)
      if opt.save_average
       % save average
       save_mgh(nanmean(cat(3,scale_metrics{:}),3),[filebase '_avg.mgh'],eye(4));
      end
      fprintf('done\n')
     end
     
     if opt.save_freesurfer && ~strcmp(opt.dim_reduction{dim_i},'none') 
      % now save as Freesurfer MGH per hemi
      fprintf('Saving as Freesurfer MGH/lh-rh...')
      % lh
      this_pos=(all_layouts{1}.hemi==1)&(all_layouts{1}.cortex<2); % choose hemi and no subcortex
      % safe
      nmri_write_mgh([filebase '_lh.mgh'],eye(4),cellfun(@(x) x(this_pos),scale_metrics,'UniformOutput',false))
      % rh
      this_pos=(all_layouts{1}.hemi==2)&(all_layouts{1}.cortex<2); % choose hemi and no subcortex
      % safe
      nmri_write_mgh([filebase '_rh.mgh'],eye(4),cellfun(@(x) x(this_pos),scale_metrics,'UniformOutput',false))
      fprintf('done\n')
     end
     
     
     if opt.save_grp_mean && ~strcmp(opt.dim_reduction{dim_i},'none') 
      % calculate a mean over all subjects (per vertex)
      fprintf('Saving group mean (suma_all)...')
      mean_map=nanmean(cat(3,scale_metrics{:}),3); % make to double (along 3rd) and mean
      % and safe
      nmri_write_mgh([filebase '_mean.mgh'],eye(4),{mean_map})
      fprintf('done\n')
     end
     
     write_log([filebase '_log'],all_subjects);
    end 
   end 
  end 
 end 
end 



%% Local functions
function metric=reduce_dim(metric,stat)
 if iscell(metric)
  metric=metric{1};
 end
 if size(metric,2)>1 && size(metric,1)>1
  % reduce
  switch(stat)
   case 'abs_mean'
    metric=abs(metric);
    metric=nanmean(metric,1);
   case 'mean'
    metric=nanmean(metric,1);
   case 'pos_mean'
    metric=metric(metric>0);
    metric=nanmean(metric,1);
   case 'neg_mean'
    metric=metric(metric<0);
    metric=nanmean(metric,1);
   case 'abs_median'
    metric=abs(metric);
    metric=median(metric,'omitnan');
   case 'median'
    metric=median(metric,'omitnan');
   case 'pos_median'
    metric=metric(metric>0);
    metric=median(metric,'omitnan');
   case 'neg_median'
    metric=metric(metric<0);
    metric=median(metric,'omitnan');
   case 'none'
    % do nothing
   otherwise
    warning(['Unkown statistical operation =' stat])
    metric=NaN(1,size(metric,1));
  end
 end
end

function write_log(fname, all_subjects)
 fid=fopen(fname,'w');
 for i=1:length(all_subjects)
  fprintf(fid,'%s\t%s\t%d\n',all_subjects{i}.id,all_subjects{i}.exam_id,i);
 end
 fclose(fid);
end

end