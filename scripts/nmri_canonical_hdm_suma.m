function [ subject ] = nmri_canonical_hdm_suma(subject, data, params)
%[ subject ] = nmri_canonical_hdm_suma(subject, data, params)
%  
% This function will use a canocical / atlas headmodel for later source
% processing,
% 
% subject   =   subject structure, see nmri_read_subject
%               could be a struct or .m file or .mat file

% written by NF 07/2021

% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct or .m/.mat file to work with')
end

% call the subject and params include
nmri_include_read_ps

% check if we have clean_dataset (maybe ICA or not)
if (isfield(params,'useICA_clean') && params.useICA_clean==1)
  if (isfield(subject,'cleanICA_dataset') &&  exist(subject.cleanICA_dataset,'file'))
   input=subject.cleanICA_dataset;
  else
   error('ICA cleaned dataset not found - Run nmri_artifactrejection first')
  end
else
 if (isfield(subject,'clean_dataset') &&  exist(subject.clean_dataset,'file'))
   input=subject.clean_dataset;
  else
   error('Cleaned dataset (w/o ICA) not found - Run nmri_artifactrejection first')
 end  
end

if (~exist('data','var') || isempty(data) )  
 disp(sprintf('Loading data = %s ...',input))
 load(input,'data');
 
 if (~exist('data','var') ) 
  error('Could not load data')
 end
end

% call central dataset selector
subject=nmri_determine_datatype(subject);

%% Get the modality-specific analysis params
[ params ] = nmri_get_modality_params( params, subject.dtype );

% set a default HDM class
if ~isfield(params,'canon_HDM')
 params.canon_HDM='nmriSUMA150_fs6'; %default
end

if isfield(subject,'active_hdm_class')
 % take from this
 hdm_class=subject.active_hdm_class;
 % take hdm type
 items=strsplit(hdm_class,'_');
 % make sure that the type is fitting
 if length(items)>1
  params.headmodel=items{end};
 end
 % fuse subject and fs class
 if length(items)>2
  putative=[items{1} '_' items{2}];
 else 
  putative=items{1};
 end
 if ~strcmp(putative,params.canon_HDM)
  hdm_class=[params.canon_HDM '_' params.headmodel];
  subject.active_hdm_class=hdm_class;
  warning(['The active HDM class set (' items{1} ') does not match the canonical HDM type set in params (' params.canon_HDM '). Switching the active class to ' hdm_class '.'])
 end
else

 hdm_class=[params.canon_HDM '_' params.headmodel];
end


%% make our output paths and dir

% we now strictly avoid mixing of different head model classes
if (~isfield(subject,'hdm_classes'))
 subject.hdm_classes=[];
end


if (~isfield(subject.hdm_classes,hdm_class))
 subject.hdm_classes.(hdm_class)=[];
end

if (~isfield(subject,'mri_dir'))
 subject.mri_dir=fullfile(subject.analysis_dir,subject.id,'mri');
end
if (~exist(subject.mri_dir,'dir'))
 mkdir(subject.mri_dir)
end

if (~isfield(subject,'proc_dir'))
 subject.proc_dir=fullfile(subject.analysis_dir,subject.id,'processed');
end
if (~exist(subject.proc_dir,'dir'))
 mkdir(subject.proc_dir)
end

% These are subject specific, so could be re-used for all exams

% commented this out, to avoid problems with subjects being moved

%if (~isfield(subject.hdm_classes.(hdm_class),'mri_ctf'))
 subject.hdm_classes.(hdm_class).mri_ctf=fullfile(subject.mri_dir,['mri_ctf_canon_' params.canon_HDM '.mat']);
%end
%if (~isfield(subject.hdm_classes.(hdm_class),'mri_seg'))
 subject.hdm_classes.(hdm_class).mri_seg=fullfile(subject.mri_dir,['mri_seg_canon_' params.canon_HDM '.mat']);
%end
%if (~isfield(subject.hdm_classes.(hdm_class),'mri_T1'))
 subject.hdm_classes.(hdm_class).mri_T1=fullfile(subject.mri_dir,['mri_T1_canon_' params.canon_HDM '.nii']);
%end
%if (~isfield(subject.hdm_classes.(hdm_class),'suma_surface'))
 subject.hdm_classes.(hdm_class).suma_surface=fullfile(subject.mri_dir,['suma_surface_canon_' params.canon_HDM '.mat']);
%end

% This is subject and exam_id specific (leadfield may differ), check if the
% currently requested headmodel is used 
%if (~isfield(subject.hdm_classes.(hdm_class),'hdm_lead')) || (isempty(regexp(subject.hdm_classes.(hdm_class).hdm_lead,['hdm_lead_canon_.*'  hdm_class '.mat$'])))
 % always set to current standard
% subject.hdm_classes.(hdm_class).hdm_lead=fullfile(subject.analysis_dir,subject.id,'processed',['hdm_lead_canon_' subject.id '_' subject.exam_id '_' hdm_class '.mat']);
%end

% new style since 02/2022 to only write subject specific lead field (and
% not a copy of the canonical headmodel)
if (~isfield(subject.hdm_classes.(hdm_class),'lead')) || (isempty(regexp(subject.hdm_classes.(hdm_class).lead,['lead_only_canon_.*'  hdm_class '.mat$'])))
 % always set to current standard
 subject.hdm_classes.(hdm_class).lead=fullfile(subject.analysis_dir,subject.id,'processed',['lead_only_canon_' subject.id '_' subject.exam_id '_' hdm_class '.mat']);
end

if (~isfield(subject.hdm_classes.(hdm_class),'hdm')) || (isempty(regexp(subject.hdm_classes.(hdm_class).hdm,['hdm_only_canon_.*'  hdm_class '.mat$'])))
 % always set to current standard
 subject.hdm_classes.(hdm_class).hdm=fullfile(subject.analysis_dir,subject.id,'processed',['hdm_only_canon_' subject.id '_' subject.exam_id '_' hdm_class '.mat']);
end

%% Now check if we know the headmodel requested

%HDMpath=[root_dir '/conf/canon_HDM/' subj_id];
% check where we have a headmodel, could be in common tool path (for large HDMs) or in the
% analyis dir

% first, check in the analysis dir
if exist(fullfile(subject.analysis_dir,'conf','canon_HDM',params.canon_HDM,['hdm_' params.headmodel '.mat']),'file')
 HDMpath=fullfile(subject.analysis_dir,'conf','canon_HDM',params.canon_HDM);
elseif exist(fullfile(getenv('NMRI_TOOLS'),'common','nmri_canonical_HDMs',params.canon_HDM,['hdm_' params.headmodel '.mat']),'file')
 HDMpath=fullfile(getenv('NMRI_TOOLS'),'common','nmri_canonical_HDMs',params.canon_HDM);
else
 error(['Could not find a canonical headmodel "' params.canon_HDM '" and type "' params.headmodel '" in the analyis dir and  common storage'])
end


% check and load HDM
hdmFile=(fullfile(HDMpath,['hdm_' params.headmodel '.mat']));

load(hdmFile); % load all that is there (could be with fem or without)

% and make symbolic link
if ~exist(subject.hdm_classes.(hdm_class).hdm,'file')
 system(['ln -s ' hdmFile ' ' subject.hdm_classes.(hdm_class).hdm] );
else
 fprintf('HDM already in place in %s\n',subject.hdm_classes.(hdm_class).hdm)
end


% check and copy suma
if ~exist(fullfile(HDMpath,['suma_all_' params.SUMA_ld '.mat']),'file')
 error(['Could not find a canonical "' params.canon_HDM '" SUMA for LD= ' params.SUMA_ld ])
end

% make symbolic link
%system(['ln -s ' fullfile(HDMpath,['suma_all_' params.SUMA_ld '.mat']) ' ' subject.hdm_classes.(hdm_class).suma_surface] );

% SUMA surface gets written always to include the canonical subject, plus
% this is small
load(fullfile(HDMpath,['suma_all_' params.SUMA_ld '.mat']),'suma_all','cfg');
canonical_subject=params.canon_HDM;
save(subject.hdm_classes.(hdm_class).suma_surface,'suma_all','canonical_subject','cfg');


% check and load MRI ctf
if ~exist(fullfile(HDMpath,'mri_ctf.mat'),'file')
 error(['Could not find a canonical "' params.canon_HDM '" MRI CTF'])
end

% make symbolic link
if ~exist(fullfile(HDMpath,'mri_ctf.mat'),'file')
 system(['ln -s ' fullfile(HDMpath,'mri_ctf.mat') ' ' subject.hdm_classes.(hdm_class).mri_ctf] );
end

% still load it for the purpose of this script
load(fullfile(HDMpath,'mri_ctf.mat'),'mri_ctf');

%save(subject.hdm_classes.(hdm_class).mri_ctf,'mri_ctf');



%% Now do the electrode alignement for EEG
if (strcmp(subject.dtype,'EEG'))
 if (~isfield(subject.hdm_classes.(hdm_class),'electrodes_aligned'))
  subject.hdm_classes.(hdm_class).electrodes_aligned=fullfile(subject.analysis_dir,subject.id,'processed',['electrodes_aligned_canon_' canonical_subject '_' subject.id '_' subject.exam_id '.mat']);
 end
 
 % determine the electrode type
 if ~isempty(regexpi(subject.detected_datatype,'EGI-.*25.?'))
  % EGI- type
  elType='egi256';
 else
  % default to 10-05
  elType='eeg1005';
 end
 
 % always check electrodes, fast and may make problems if channels were
 % changed
 
 
 %if (exist(subject.hdm_classes.(hdm_class).electrodes_aligned,'file'))
 % disp('Found aligend EEG electrodes - loading')
 % load(subject.hdm_classes.(hdm_class).electrodes_aligned,'elec_aligned','elec_present','elec_missing')
 %else
  if (~isfield(subject,'elec_file') || ~isfield(subject,'elec_fiducials') || isempty(subject.elec_file))
   error('Electrodes information missing in subject struct, should be auto-determined by dataset selector. Check dataset-mappings files.')
  end
  
  
  % check and load electrodes



  if ~exist(fullfile(HDMpath,['electrodes_aligned_' elType '.mat']),'file')
   error(['Could not find a canonical "' params.canon_HDM '" electrodes for type=' elType ])
  end
  load(fullfile(HDMpath,['electrodes_aligned_' elType '.mat']),'elec_aligned');


  % sort electrodes into missing and present
  elec_present=[];
  elec_present.unit=elec_aligned.unit;
  elec_present.label={};
  elec_present.chantype={};
  elec_present.chanunit={};
  elec_present.chanpos=[];
  elec_present.elecpos=[];
  elec_missing=[];
  elec_missing.unit=elec_aligned.unit;
  elec_missing.label={};
  elec_missing.chantype={};
  elec_missing.chanunit={};
  elec_missing.chanpos=[]; 
  elec_missing.elecpos=[]; 
  % make UPPER case labels -- seems important
  elec_aligned.label=upper(elec_aligned.label);
  
  % deal with T3/T4, T5/T6 10-20 mess...
  if strcmp(elType,'eeg1005')
   if any(strcmpi(data.label,'T3')) && any(strcmpi(elec_aligned.label,'T7'))
    elec_aligned.label{strcmpi(elec_aligned.label,'T7')}='T3';
   end
   if any(strcmpi(data.label,'T4')) && any(strcmpi(elec_aligned.label,'T8'))
    elec_aligned.label{strcmpi(elec_aligned.label,'T8')}='T4';
   end
   if any(strcmpi(data.label,'T5')) && any(strcmpi(elec_aligned.label,'P7'))
    elec_aligned.label{strcmpi(elec_aligned.label,'P7')}='T5';
   end
   if any(strcmpi(data.label,'T6')) && any(strcmpi(elec_aligned.label,'P8'))
    elec_aligned.label{strcmpi(elec_aligned.label,'P8')}='T6';
   end
  end 
  
  
  for i=1:length(elec_aligned.label)
   match_elec=find(strcmp(data.label,elec_aligned.label{i})); % NOT do case insens.
   
   % deal with channels marked as BAD, as of 10/2019
   if isfield(data,'bad_channels') && any(strcmp(elec_aligned.label{i},data.bad_channels))
    % electrode is present, but bad  -treat as missing
    match_elec=[];
   end
   
   % check if already present, could happen with 1005
   if any(strcmp(elec_present.label,elec_aligned.label{i})) || any(strcmp(elec_missing.label,elec_aligned.label{i}))
    % do not add to either
    continue
   end
   
   if isempty(match_elec)
    elec_missing.label{end+1,1}=elec_aligned.label{i};
    elec_missing.chantype{end+1,1}=elec_aligned.chantype{i,1};
    elec_missing.chanunit{end+1,1}=elec_aligned.chanunit{i,1};
    elec_missing.chanpos(end+1,1:3)=elec_aligned.chanpos(i,1:3); 
    elec_missing.elecpos(end+1,1:3)=elec_aligned.elecpos(i,1:3); 
   else
    elec_present.label{end+1,1}=elec_aligned.label{i,1};
    elec_present.chantype{end+1,1}=elec_aligned.chantype{i,1};
    elec_present.chanunit{end+1,1}=elec_aligned.chanunit{i,1};
    elec_present.chanpos(end+1,1:3)=elec_aligned.chanpos(i,1:3); 
    elec_present.elecpos(end+1,1:3)=elec_aligned.elecpos(i,1:3);
   end
  end
  
  canonical_subject=params.canon_HDM;
  save(subject.hdm_classes.(hdm_class).electrodes_aligned,'elec_aligned','elec_present','elec_missing','canonical_subject');
 %end
end
%% Make sure we have OpenMEEG binaries, if needed
switch params.headmodel
 case 'openmeeg'
   % OpenMEEG BEM head model
   
   % this is an external toolbox that needs to be in the path to work
   [ret,~]=system('om_assemble');
   if ret==0
    fprintf('Found OpenMEEG\n')
   else
    % for GWDG/Göttingen try to fix it
    setenv('PATH', [ getenv('NMRI_TOOLS')  '/OpenMEEG/conda/bin:' getenv('PATH')]);
    setenv('LD_LIBRARY_PATH', [ getenv('NMRI_TOOLS')  '/OpenMEEG/conda/lib:' getenv('LD_LIBRARY_PATH')]);
    % and try again
    [ret,~]=system('om_assemble');
    if ret==0
     fprintf('Found OpenMEEG after adding it to the PATH\n')
    else
     error('OpenMEEG binaries not found. Make sure these are in the PATH.')
    end
   end
      
end

%% now make leadfield - this should be done always, since we may have different sensors
cfg = [];
if (strcmp(subject.dtype,'MEG') && isfield(data,'grad'))
 cfg.grad    = data.grad;                  % sensor information
elseif (strcmp(subject.dtype,'EEG') && exist('elec_aligned','var'))
 

 cfg.elec    = elec_present;     % electrode positions
else
 error('No gradiometers or electrodes defined')
end
cfg.channel = data.label;                 % the used channels

cfg.sourcemodel.pos = double(suma_all.pos);              % source points, some Fieldtrip functions need this to be double
cfg.headmodel = hdm; 
cfg.normalize = 'no';                    % set to no as of 08/March/2017, recommended by Fieldtrip mailing list
 
[leadfield cfg]= ft_prepare_leadfield(cfg);


%% save leadfield incl. subject info
canonical_subject=params.canon_HDM;

if exist('fem','var')
 save(subject.hdm_classes.(hdm_class).lead,'leadfield','bnd','subject','fem','canonical_subject','-v7.3'); % changed to v7.3 to cover larger files
else
 save(subject.hdm_classes.(hdm_class).lead,'leadfield','bnd','subject','canonical_subject','-v7.3'); % changed to v7.3 to cover larger files
end

%% if we got this far, place a stamp for completed hdm-SUMA
subject=nmri_stamp_subject(subject,['hdmSUMA_' canonical_subject '_' params.headmodel],params);
 

%% visually check the alignement
disp('Generating brain / sensors figure...')
hFig=figure('Position',[0,0,800,800],'Visible','off');
hFig=figure('Position',[0,0,800,800],'Visible','on');
hold on     % plot all objects in one figure
switch params.headmodel
 % note, if complied the eval include of color definition does not work,
 % set fixed colors now
 case 'singleshell'
  % usual MEG model
  ft_plot_headmodel(hdm, 'facecolor',[0.781 0.762 0.664],'edgecolor', 'none','facealpha',0.6);
  ft_plot_mesh(bnd(3),'facecolor',[1 0.878 0.741],'edgecolor','none','facealpha',0.2);  
 case 'dipoli'
  ft_plot_mesh(bnd(1),'facecolor',[0.781 0.762 0.664],'edgecolor','none','facealpha',0.6);  
  ft_plot_mesh(bnd(2),'facecolor',[1 1 1],'edgecolor','none','facealpha',0.2);  
  ft_plot_mesh(bnd(3),'facecolor',[1 0.878 0.741],'edgecolor','none','facealpha',0.2);
 case 'openmeeg'
  ft_plot_mesh(bnd(1),'facecolor',[0.781 0.762 0.664],'edgecolor','none','facealpha',0.6);  
  ft_plot_mesh(bnd(2),'facecolor',[1 1 1],'edgecolor','none','facealpha',0.2);  
  ft_plot_mesh(bnd(3),'facecolor',[1 0.878 0.741],'edgecolor','none','facealpha',0.2);  
 case 'simbio'
  ft_plot_mesh(fem, 'surfaceonly', 'yes','facecolor',[0.781 0.762 0.664],'edgecolor', 'none','facealpha',0.3);
end
% plot the leadfield points
% note, for compiled run, needs facecolor to be set 
ft_plot_mesh(leadfield.pos(leadfield.inside,:),'vertexcolor','b','facecolor','none');
% show not "inside" in large and red
ft_plot_mesh(leadfield.pos(~leadfield.inside,:),'vertexcolor','r','vertexsize',15,'facecolor','none');

% plot fiducials 
nas=[mri_ctf.cfg.fiducial.nas 1]*mri_ctf.transform';
lpa=[mri_ctf.cfg.fiducial.lpa 1]*mri_ctf.transform';
rpa=[mri_ctf.cfg.fiducial.rpa 1]*mri_ctf.transform';
ft_plot_mesh(nas(1:3),'vertexcolor',[1 1 1],'vertexsize',30,'vertexmarker','x','facecolor','none');
ft_plot_mesh(lpa(1:3),'vertexcolor',[0 1 0],'vertexsize',30,'vertexmarker','x','facecolor','none');
ft_plot_mesh(rpa(1:3),'vertexcolor',[1 0 0],'vertexsize',30,'vertexmarker','x','facecolor','none');

if (strcmp(subject.dtype,'MEG') && isfield(data,'grad'))
 ft_plot_sens(data.grad,'style', 'g*');
 clipping=[500 700 150 0];
elseif (strcmp(subject.dtype,'EEG') && exist('elec_aligned','var'))
 % plot present always
 ft_plot_mesh(elec_present.chanpos,'vertexmarker','k.','vertexcolor',[0 0 0],'vertexsize',12,'facecolor','none');
 % missing needs to be dealt with seperately for 10-05 (many channels
 % usually not there)
 if strcmp(elType,'eeg1005')
  % show only those as missing that are bad channels
  if isfield(data,'bad_channels')
   badIdx=zeros(length(elec_missing.label),1);
   for ii=1:length(data.bad_channels)
    sel=find(strcmp(elec_missing.label,data.bad_channels{ii}));
    badIdx(sel)=1;
   end
   if any(badIdx)
    ft_plot_mesh(elec_missing.chanpos(badIdx,:),'vertexmarker','k.','vertexcolor',[0.8 0.5 0.5],'vertexsize',15,'facecolor','none');
   end
  end
 else
  % default to show all missing
  ft_plot_mesh(elec_missing.chanpos,'vertexmarker','k.','vertexcolor',[0.8 0.5 0.5],'vertexsize',15,'facecolor','none');
 end
 clipping=[800 800 0 0];
end
hold off
camlight('headlight');

 
% Export figure if wanted
if (isfield(params,'QC_plots') && params.QC_plots==1)
 % make the QC dir
 if isfield(subject,'QCdir')
  QCdir=subject.QCdir;
 else
  QCdir=fullfile(subject.analysis_dir,subject.id,'QC');
 end
 if (~exist(QCdir,'dir'))
  mkdir(QCdir)
 end
 
 hFig=nmri_rotate_print(hFig,[0 0;180 0;90 0;0 90],fullfile(QCdir,['hdm_sensors_plot_canon_' subject.id '_' subject.exam_id '_' hdm_class]),[subject.id ' - ' params.canon_HDM],clipping);
end
 
% only show after save
if isdeployed || ~usejava('desktop')
 % headless, do not show
 close(hFig)
else
 title(subject.id);
 view([0 0])
 set(hFig,'Visible','on')
end
 
disp('...canonical headmodelling done')
disp('Check the alignement now and proceed to nmri_processing when happy'); 
 
end

