function [ subject ] = nmri_read_markers( subject, params )
%[ subject ] = nmri_read_markers( subject, params )
%   Will check spike/event markers and put into subject struct

% call the subject and params include
nmri_include_read_ps

%% Get the modality-specific analysis params
if isfield(subject,'dataset_mapping')
 [ params ] = nmri_get_dataset_params( params, subject.dataset_mapping );
end
[ params ] = nmri_get_modality_params( params, subject.dtype );

%% Safe subject in
subject_in=subject;

%% check the cleaned version
if (exist(subject.dws_filt_dataset,'file'))
 load(subject.dws_filt_dataset,'data');   
else
 error('We need the downsampled / filtered dataset to run marker reading. Use nmri_preproc first')
end

%% use the new more flexible event / markings reader
allPatterns={};
allFiles={};

if isfield(params,'readEvents')

 
 subject.evt_markerFile={};
 subject.evt_timings_sample=[];
 subject.evt_timings_seconds=[];
 subject.evt_IDs={};

 subject.info_markerFile={};
 subject.info_timings_sample=[];
 subject.info_timings_seconds=[];
 subject.info_IDs={};


 % now loop for all classes of events
 allEvts=fieldnames(params.readEvents);
 for evC=1:length(allEvts)

  % set path to look 
  pattern=params.readEvents.(allEvts{evC}).Path;

  % translate our placeholders
  pattern=strrep(pattern,'<subject_id>',subject.id);
  pattern=strrep(pattern,'<exam_id>',subject.exam_id);
  pattern=strrep(pattern,'<raw_dataset>',subject.raw_dataset); % dataset, e.g. of EDF reader
  pattern=strrep(pattern,'<raw_datadir>',fileparts(subject.raw_dataset)); % path of raw dataset
  pattern=strrep(pattern,'<analysis_dir>',subject.analysis_dir);
 
  % make cell array 

  if ischar(pattern)
   pattern={pattern};
  end
  
  allPatterns=[allPatterns pattern];
  
   
  % loop over paths
  found_files={};
  for iP=1:length(pattern)
   eventFileThis=dir(pattern{iP});
   for iF=1:length(eventFileThis)
    if ~isfield(params.readEvents.(allEvts{evC}),'regexp') || ~isempty(regexp(eventFileThis(iF).name,params.readEvents.(allEvts{evC}).regexp))
     found_files{end+1,1}=fullfile(eventFileThis(iF).folder,eventFileThis(iF).name);
    end
   end
  end
  
  if (~isempty(found_files))
   % found at least one
   if (length(found_files)>1)
    % multiple, so choose
    manselect=listdlg('Name','Choose Event Marker File','PromptString',{'Found >1 possibility, please select (can be multiple)',['Subject=' subject.id  ', ExamID=' subject.exam_id ', EventType=' allEvts{evC}]},'SelectionMode','multiple','ListSize',[700 (50+(length(found_files)*10))],'ListString',found_files);
    if (~isempty(manselect))
     found_files=found_files(manselect);
    else
     error('No selection was made -- stopping here')
    end
   end
   
   
   % now loop over found files
   for ff=1:length(found_files)
    fprintf('Now parsing file=%s for events\n',found_files{ff})
    if isfield(params.readEvents.(allEvts{evC}),'EventsTextscan')
     % Do textscan, usually fine for text files
     fprintf('Using textscan parsing\n')
     fid = fopen(found_files{ff}); % file ID
     eval(params.readEvents.(allEvts{evC}).EventsTextscan);
     fclose(fid);
     infoIDs=[];
     infoTimes=[];
     % all are spikes
    elseif isfield(params.readEvents.(allEvts{evC}),'EventsTableRead')
     % Do table import, flexible for many formats (CSV, xls,...)
     fprintf('Using readtable parsing\n')
     read_table = readtable(found_files{ff});
     eval(params.readEvents.(allEvts{evC}).EventsTableRead); 
     % try to sort
     if ~isfield(params.readEvents.(allEvts{evC}),'Class') || strcmp(params.readEvents.(allEvts{evC}).Class,'Auto')
      % auto-determine spikes
      isSpike=~cellfun(@isempty,regexpi(spikeIDs,'spk')) | ~cellfun(@isempty,regexpi(spikeIDs,'spike'));
      infoTimes=spikeTimes(~isSpike);
      infoIDs=spikeIDs(~isSpike);
      spikeTimes(~isSpike)=[];
      spikeIDs(~isSpike)=[];
     elseif strcmp(params.readEvents.(allEvts{evC}).Class,'Spikes')
      infoIDs=[];
      infoTimes=[];
      % all are spikes
     elseif strcmp(params.readEvents.(allEvts{evC}).Class,'Infos')
      infoIDs=spikeIDs;
      infoTimes=spikeTimes;
      spikeIDs=[];
      spikeTimes=[];
      % all are infos
     end
    elseif isfield(params.readEvents.(allEvts{evC}),'EventsFunction')
     % call a specific function, e.g. EGI xml oder EDF
     fprintf('Using function (%s) parsing\n',params.readEvents.(allEvts{evC}).EventsFunction)
     [ this_events, this_infos ]=feval(params.readEvents.(allEvts{evC}).EventsFunction,found_files{ff},subject);
     spikeTimes=cell2mat(this_events.evt_timings_seconds);
     spikeIDs=this_events.evt_IDs;
     infoTimes=cell2mat(this_infos.evt_timings_seconds);
     infoIDs=this_infos.evt_IDs;
    else
     error(['No extraction function specified for event class=' allEvts{evC} '. Check analysis_params.m.' ])
    end
   

    spikeTimes = double(spikeTimes);
    nSpikes = length(spikeTimes);
    % check for spikeNames
    if ~exist('spikeIDs','var') || length(spikeIDs)~=nSpikes
     spikeIDs=repmat({'Spike'},nSpikes,1);
    end
    
    infoTimes = double(infoTimes);
    nInfos = length(infoTimes);
    if ~exist('infoIDs','var') || length(infoIDs)~=nInfos
     infoIDs=repmat({'Marking'},nInfos,1);
    end

    % convert the time into seconds of the spikes markers
    spikeTimesSeconds = spikeTimes/params.readEvents.(allEvts{evC}).EventsTimebase;
    infoTimesSeconds = infoTimes/params.readEvents.(allEvts{evC}).EventsTimebase;

    % sample number that has the spike markings
    spikes = round(spikeTimesSeconds*data.fsample)+1;
    infos = round(infoTimesSeconds*data.fsample)+1;

    subject.evt_markerFile=[ subject.evt_markerFile found_files{ff}];
    subject.evt_timings_sample=[subject.evt_timings_sample spikes];
    subject.evt_timings_seconds=[subject.evt_timings_seconds spikeTimesSeconds];
    subject.evt_IDs=[subject.evt_IDs spikeIDs];

    subject.info_markerFile=[ subject.info_markerFile found_files{ff} ];
    subject.info_timings_sample=[ subject.info_timings_sample infos];
    subject.info_timings_seconds=[ subject.info_timings_seconds infoTimesSeconds];
    subject.info_IDs=[subject.info_IDs infoIDs];
  
    % remove old not found marker
    if isfield(subject,'evt_markerFile_notFound')
     subject=rmfield(subject,'evt_markerFile_notFound');
    end
    
   end
   
   allFiles=[allFiles found_files];
   
   fprintf('Found %d Events and %d Info Annotations\n',length(spikeIDs),length(infoIDs))
  end
 end
else
 disp('No event path set, not parsing the file / marker')
end 

%% Now check if we expect to have found something

if isfield(params,'rejectEvents') && params.rejectEvents==1
 if ~isfield(subject,'evt_markerFile') || isempty(subject.evt_markerFile)  || isempty(subject.evt_timings_sample)
  % no marker found, give warning if not a control
  if (~strcmp(subject.id(1:2),'C.'))
   h=msgbox({'No spike markings found: please make sure that this patient has no spikes!','',['Search pattern:' strjoin(allPatterns,' ')],'','Terminate processing, if needed'},'No events found');
   uiwait(h)
   subject.evt_markerFile_notFound=pattern;
  else
   disp('This seems to be a control subject and no event marks were found - this is probably okay')
  end
 end
end


%% Parse stimuli if needed
dws_saved=0;
if isfield(params,'parse_stimuli') && params.parse_stimuli==1
 fprintf('Parseing stimuli...\n')
 data_out=nmri_parse_stimuli(subject, data, params);
 
 % and write out to dataset, if needed
 if ~isfield(data,'trial_markings') || ~isequal(data_out.trial_markings,data.trial_markings)
  fprintf('Saving new dws_filt dataset with stimuli info...\n')
  data=data_out;
  save(subject.dws_filt_dataset,'data','subject');
  dws_saved=1;
 end
 
 
 % and update the other datasets as well
 if isfield(subject,'cleanICA_dataset') && exist(subject.cleanICA_dataset,'file')
  data_orig=load(subject.cleanICA_dataset,'data','subject');
  if ~isfield(data_orig.data,'trial_markings') || size(data_orig.data.trial_markings,2)<4 ||...
    ~isequal(data_orig.data.trial_markings(:,4),data_out.trial_markings(:,4))
   qcfg=[];
   qcfg.question={'Changes in the provocations markings have been found.','',...
   'Do you want to overwrite the markings in CleanedICA?',...
   'Note, this will erase manually placed provocations markings.'};
   qcfg.title='Overwrite?';
   qcfg.options={'Yes','No'};
   qcfg.mode='buttons';
   response=nf_uidialog(qcfg);
   if strcmpi(response,'Yes')
    fprintf('Saving updated ICA cleaned dataset with stimuli info...\n')
    data=data_orig.data;
    % overwrite stimuli
    data.trial_markings(:,4)=data_out.trial_markings(:,4);
    save(subject.cleanICA_dataset,'data','-append');
   end
  end
  clear data_orig
 end
 
 if isfield(subject,'clean_dataset') && exist(subject.clean_dataset,'file')
  data_orig=load(subject.clean_dataset,'data','subject');
  if ~isfield(data_orig.data,'trial_markings') || size(data_orig.data.trial_markings,2)<4 ||...
    ~isequal(data_orig.data.trial_markings(:,4),data_out.trial_markings(:,4))
   qcfg=[];
   qcfg.question={'Changes in the provocations markings have been found.','',...
   'Do you want to overwrite the markings in CleanedICA?',...
   'Note, this will erase manually placed provocations markings.'};
   qcfg.title='Overwrite?';
   qcfg.options={'Yes','No'};
   qcfg.mode='buttons';
   response=nf_uidialog(qcfg);
   if strcmpi(response,'Yes')
    fprintf('Saving updated cleaned dataset with stimuli info...\n')
    data=data_orig.data;
    % overwrite stimuli
    data.trial_markings(:,4)=data_out.trial_markings(:,4);
    save(subject.clean_dataset,'data','subject');
   end
  end
  clear data_orig
 end
end


%% Now write out things
% save original subject
subject_local=subject;
if (~isequal(subject_in, subject))
 % save again, to be keep this info in mat file, unless done already
 if ~dws_saved
  save([subject.dws_filt_dataset],'subject','-append');
 end
 
 % add to other datasets, just in case
 if ~isfield(subject,'evt_markerFile_notFound') 
  % start with most advanced
  if isfield(subject,'hdm_lead') && exist(subject.hdm_lead,'file')
   load(subject.hdm_lead,'subject')
   if ~isfield(subject,'evt_timings_seconds') || ~isequal(subject.evt_timings_seconds,spikeTimesSeconds) || ~isfield(subject,'info_timings_seconds') || ~isequal(subject.info_timings_seconds,infoTimesSeconds) ...
     || ~isfield(subject,'evt_IDs') || ~isequal(subject.evt_IDs,spikeIDs) || ~isfield(subject,'info_IDs')  || ~isequal( subject.info_IDs,infoIDs)
    % any change, so update
    subject.evt_markerFile=allFiles;
    subject.evt_timings_sample=spikes;
    subject.evt_timings_seconds=spikeTimesSeconds;
    subject.evt_IDs=spikeIDs;
    subject.info_markerFile=allFiles;
    subject.info_timings_sample=infos;
    subject.info_timings_seconds=infoTimesSeconds;
    subject.info_IDs=infoIDs;
    fprintf('Updating file = %s\n',subject.hdm_lead)
    save([subject.hdm_lead],'subject','-append');
   end
  end
  
  if isfield(subject,'cleanICA_dataset') && exist(subject.cleanICA_dataset,'file')
   load(subject.cleanICA_dataset,'subject')
   if ~isfield(subject,'evt_timings_seconds') || ~isequal(subject.evt_timings_seconds,spikeTimesSeconds) || ~isfield(subject,'info_timings_seconds') || ~isequal(subject.info_timings_seconds,infoTimesSeconds) ...
     || ~isfield(subject,'evt_IDs') || ~isequal(subject.evt_IDs,spikeIDs) || ~isfield(subject,'info_IDs')  || ~isequal( subject.info_IDs,infoIDs)
     qcfg=[];
     qcfg.question={'Changes in the events / spikes markings have been found.','',...
    'Do you want to overwrite the markings in CleanedICA?',...
    'Note, this will erase manually placed spike / event markings.'};
    qcfg.title='Overwrite?';
    qcfg.options={'Yes','No'};
    qcfg.mode='buttons';
    response=nf_uidialog(qcfg);
    if strcmpi(response,'Yes')
     % any change, so update
     subject.evt_markerFile=allFiles;
     subject.evt_timings_sample=spikes;
     subject.evt_timings_seconds=spikeTimesSeconds;
     subject.evt_IDs=spikeIDs;
     subject.info_markerFile=allFiles;
     subject.info_timings_sample=infos;
     subject.info_timings_seconds=infoTimesSeconds;
     subject.info_IDs=infoIDs;
     fprintf('Updating file = %s\n',subject.cleanICA_dataset)
     save([subject.cleanICA_dataset],'subject','-append');
    end
   end
  end
  
  if isfield(subject,'ICA_components') && exist(subject.ICA_components,'file')
   load(subject.ICA_components,'subject')    
   if ~isfield(subject,'evt_timings_seconds') || ~isequal(subject.evt_timings_seconds,spikeTimesSeconds) || ~isfield(subject,'info_timings_seconds') || ~isequal(subject.info_timings_seconds,infoTimesSeconds) ...
     || ~isfield(subject,'evt_IDs') || ~isequal(subject.evt_IDs,spikeIDs) || ~isfield(subject,'info_IDs')  || ~isequal( subject.info_IDs,infoIDs)
    % any change, so update
    subject.evt_markerFile=allFiles;
    subject.evt_timings_sample=spikes;
    subject.evt_timings_seconds=spikeTimesSeconds;
    subject.evt_IDs=spikeIDs;
    subject.info_markerFile=allFiles;
    subject.info_timings_sample=infos;
    subject.info_timings_seconds=infoTimesSeconds;
    subject.info_IDs=infoIDs;
    fprintf('Updating file = %s\n',subject.ICA_components)
    save([subject.ICA_components],'subject','-append');
   end
  end
  
  
  if isfield(subject,'clean_dataset') && exist(subject.clean_dataset,'file')
   load(subject.clean_dataset,'subject')    
   if ~isfield(subject,'evt_timings_seconds') || ~isequal(subject.evt_timings_seconds,spikeTimesSeconds) || ~isfield(subject,'info_timings_seconds') || ~isequal(subject.info_timings_seconds,infoTimesSeconds) ...
     || ~isfield(subject,'evt_IDs') || ~isequal(subject.evt_IDs,spikeIDs) || ~isfield(subject,'info_IDs')  || ~isequal( subject.info_IDs,infoIDs)
    qcfg=[];
    qcfg.question={'Changes in the events / spikes markings have been found.','',...
    'Do you want to overwrite the markings in the cleaned dataset?',...
    'Note, this will erase manually placed spike / event markings.'};
    qcfg.title='Overwrite?';
    qcfg.options={'Yes','No'};
    qcfg.mode='buttons';
    response=nf_uidialog(qcfg);
    if strcmpi(response,'Yes')
     % any change, so update
     subject.evt_markerFile=allFiles;
     subject.evt_timings_sample=spikes;
     subject.evt_timings_seconds=spikeTimesSeconds;
     subject.evt_IDs=spikeIDs;
     subject.info_markerFile=allFiles;
     subject.info_timings_sample=infos;
     subject.info_timings_seconds=infoTimesSeconds;
     subject.info_IDs=infoIDs;
     fprintf('Updating file = %s\n',subject.clean_dataset)
     save([subject.clean_dataset],'subject','-append');
    end
   end
  end


 else
  % not file found -- still save
  if isfield(subject,'hdm_lead') && exist(subject.hdm_lead,'file')
   load(subject.hdm_lead,'subject')
   if ~isfield(subject,'evt_markerFile_notFound') || ~isequal(subject.evt_markerFile_notFound,pattern)
    % any change, so update
    subject.evt_markerFile_notFound=pattern;
    fprintf('Updating file = %s\n',subject.hdm_lead)
    save([subject.hdm_lead],'subject','-append');
   end
  end 
  
  if isfield(subject,'cleanICA_dataset') && exist(subject.cleanICA_dataset,'file')
   load(subject.cleanICA_dataset,'subject')
   if ~isfield(subject,'evt_markerFile_notFound') || ~isequal(subject.evt_markerFile_notFound,pattern)
    % any change, so update
    subject.evt_markerFile_notFound=pattern;
    fprintf('Updating file = %s\n',subject.cleanICA_dataset)
    save([subject.cleanICA_dataset],'subject','-append');
   end
  end 
  
  if isfield(subject,'clean_dataset') && exist(subject.clean_dataset,'file')
   load(subject.clean_dataset,'subject')    
   if ~isfield(subject,'evt_markerFile_notFound') || ~isequal(subject.evt_markerFile_notFound,pattern)
    % any change, so update
    subject.evt_markerFile_notFound=pattern;
    fprintf('Updating file = %s\n',subject.clean_dataset)
    save([subject.clean_dataset],'subject','-append');
   end
  end
 end
end

% get back original subject
subject=subject_local;
  
% if we got this far, place a stamp for completed event reading
subject=nmri_stamp_subject(subject,'readingevents',params);
 


end

