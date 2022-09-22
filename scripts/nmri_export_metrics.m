function [ all_subjects ] = nmri_export_metrics( all_subjects, opt, params )
%[ all_subjects ] = nmri_export_metrics( all_subjects, opt, params )
% Will save any number of metrics / frequency bands (after scaling / 
% dim_reduction/ abs) in a 4D .mgh file. The order will be as in the 
% all_subject cell-array of structs
%

% opt               = struct array of options
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
%  .save_freesurfer = true / false, conform to freesurfer convention
%                     i.e. lh/rh and only cortex, default: false
%  .save_fsaverage  = true / false, write resampled to fsaverage, also
%                     lh/rh and only cortex, default: false
%  .save_average    = true / false, write average (for all specified 
%                     save options), default: false
%  .save_mgh        = save as .mgh files - for Freesurfer / PALM
%                     default: true
%  .save_mat        = save as .mat files - for Matlab 
%                     default: true
%  .save_global     = save the global value (1 value/per subject) as .mat file 
%                     default: true
%  .save_grp_mean   = save the mean per map/variation as .mgh
%                     default: false
%  .output          = output dir, will place 4D .mgh files there 
%                     default: <analysis_dir>/export/<date>
%  .hdm_class       = optional, specify a certain type of headmodel to use.
%                     If not set, determine from lastes used/active_hdm in
%                     subject struct.
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


if ~isfield(opt,'save_fsaverage') || ~islogical(opt.save_fsaverage ) || isempty(opt.save_fsaverage )
 opt.save_fsaverage=false;
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
if ~isfield(opt,'save_global') || ~islogical(opt.save_global ) || isempty(opt.save_global )
 opt.save_global=true;
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

suma_ld='';

% we may have a special HDM class to us
if isfield(opt,'hdm_class') && ~isempty(opt.hdm_class)
 active_hdm_class=opt.hdm_class;
else
 active_hdm_class='';
end
% check that all subjects have the same modality, headmodel and SUMA
for n=1:length(all_subjects)
 % process active subject
 subject=all_subjects{n}; 
 % get modality params
 [ this_params ] = nmri_get_modality_params( params, subject.dtype );
 
 % get the hdm_class, usually from 1st subject
 if isempty(active_hdm_class)
  if isfield(subject,'active_hdm_class')
   active_hdm_class=subject.active_hdm_class;
  elseif isfield(params,'active_hdm_class')
   active_hdm_class=params.active_hdm_class;
  else
   active_hdm_class=['individual_',this_params.headmodel];
  end 
 end
 

 % do some checks
 
 % check if we have processing done of this class
 if isfield(subject,'hdm_classes') && isfield(subject.hdm_classes,active_hdm_class) && isfield(subject.hdm_classes.(active_hdm_class),'stats')
  % then use these files from now on
  subject.stats=subject.hdm_classes.(active_hdm_class).stats;
  subject.suma_surface=subject.hdm_classes.(active_hdm_class).suma_surface;
  % fix the HDM type
  all_subjects{n}.active_hdm_class=active_hdm_class;
 else
  % old style 
  % make empty struct
  subject.hdm_classes.(active_hdm_class)=[];
 end
 
 if ~isfield(subject.hdm_classes.(active_hdm_class),'stats') || ~exist(subject.hdm_classes.(active_hdm_class).stats,'file') || ~isfield(subject.stamps,['processing_' active_hdm_class])
  error(['Subject ID=' subject.id ' does not seem to be processed for HDM= ' active_hdm_class ' - cannot continue. Process this subject or remove from subject list.'])
 end
  
 % check SUMA template for this one
 if isfield(subject.stamps.(['processing_' active_hdm_class]),'params') && isfield(subject.stamps.(['processing_' active_hdm_class]).params,'SUMA_ld')
  this_suma_ld=subject.stamps.(['processing_' active_hdm_class]).params.SUMA_ld;
  if isempty(suma_ld)
   suma_ld=this_suma_ld;
  end
  if ~strcmp(this_suma_ld,suma_ld)
   error(['Mismatch of SUMA-ld for subject=' subject.id])
  end
 else
  % we do not seem to have a param set in stamp...this shoudl actually not
  % happen, but seems to be the case sometimes
  fprintf('No params found in stamp of %s, %s\n',subject.id, subject.exam_id)
 end
end

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
  [all_metrics, all_suma, subjects_back, back_active_hdm_class]=nmri_load_metrics(all_subjects,opt.metrics{m},freq_n(f),params);
  
  if ~isequal(subjects_back,all_subjects)
   warning(['Could not find the metric ' opt.metrics{m} '/' freq_txt{f} ' for all subjects -- skipping'])
   continue
  end
  
  if ~strcmp(back_active_hdm_class,active_hdm_class)
   error(['Mismatch of HDM class of metric loader (' back_active_hdm_class ') and wrapper scrip (' active_hdm_class '). Investigate.'])
  end
  
  % get all suma msk, start with 1st subject
  all_msk=all_suma{1}.msk;
  for si=2:length(all_suma)
   if isequal(size(all_msk),size(all_suma{si}.msk))
    all_msk=all_msk&all_suma{si}.msk;
   else
    error(['Problem with the SUMA vertex count for ID=' all_subjects{si}.id])
   end
  end
  % remove from all
  for mi=1:length(all_metrics)
   if size(all_metrics{mi},2)>1 && size(all_metrics{mi},1)>1
    % two dims
    all_metrics{mi}(~all_msk,~all_msk)=NaN;  
   else
    % one dim
    all_metrics{mi}(~all_msk)=NaN;   
   end
  end

  % save mask
  if opt.save_mgh
   save_mgh(all_msk,fullfile(outputdir,'all_suma_msk.mgh'),eye(4));
  end
  if opt.save_mat
   save(fullfile(outputdir,'all_suma_msk.mat'),'all_msk','all_subjects');
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

     filebase=fullfile(outputdir,[ opt.metrics{m} '_' freq_txt{f} '_' global_txt '_' scale_txt dim_txt '_' active_hdm_class '_N' num2str(length(all_subjects))]);
     if opt.save_mat
      fprintf('Saving as .mat...')
      save([filebase '.mat'],'scale_metrics','all_subjects','active_hdm_class')
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
      this_pos=(all_suma{1}.hemi==1)&(all_suma{1}.cortex<2); % choose hemi and no subcortex
      % safe
      nmri_write_mgh([filebase '_lh.mgh'],eye(4),cellfun(@(x) x(this_pos),scale_metrics,'UniformOutput',false))
      % rh
      this_pos=(all_suma{1}.hemi==2)&(all_suma{1}.cortex<2); % choose hemi and no subcortex
      % safe
      nmri_write_mgh([filebase '_rh.mgh'],eye(4),cellfun(@(x) x(this_pos),scale_metrics,'UniformOutput',false))
      fprintf('done\n')
     end
     
     if opt.save_fsaverage && ~strcmp(opt.dim_reduction{dim_i},'none') 
      % map to fsaverage
      fprintf('Map to fsaverage...')
      [ scale_metrics_fs_lh, scale_metrics_fs_rh ] = nmri_map2fsaverage( all_subjects, scale_metrics );
      % and safe
      nmri_write_mgh([filebase '_fsaverage_lh.mgh'],eye(4),scale_metrics_fs_lh)
      nmri_write_mgh([filebase '_fsaverage_rh.mgh'],eye(4),scale_metrics_fs_rh) 
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
     
     if opt.save_global 
      % calculate a mean global over all subjects (per vertex)
      fprintf('Saving global values (reduced by %s) as .csv...',opt.dim_reduction{dim_i})
      switch(opt.dim_reduction{dim_i})
       case 'mean'
        global_metrics=cellfun(@nanmean,scale_metrics);
       case 'median'
        global_metrics=cellfun(@(x) median(x,'omitnan'),scale_metrics);
       case 'rms'
        global_metrics=cellfun(@rms,scale_metrics);
       otherwise
        global_metrics=[];
      end
    
      csvwrite([filebase '_global' dim_txt '.csv'],global_metrics)
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