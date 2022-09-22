function  [ events, infos ] = nmri_read_event_edf( edf_file, subject )
% [ events, infos ] = nmri_read_event_edf( subject )
%   Reads in EDF annotations from 


if ~exist('subject','var') || ~isfield(subject,'hdr') || ~isfield(subject.hdr,'orig') || ~isfield(subject.hdr.orig,'annotation')
 error ('Subject not set or HDR not found or not EDF(+) with annotations')
end

if ~exist(edf_file,'file')
 error(['EDF file ' subject.raw_dataset ' not found'])
end

all_evts=ft_read_event(edf_file,'header',subject.hdr);

% need a special tweak for dealing with discontinous EDF+ files
% at least NATUS seems to write out an empty annotation for each second
% hence we can calculate a corrected time by adding these (since Fieldtrip
% reads data as continuous
% checked also for Nicolet/SPZ


% now read (spike) events (SPK1)
events.evt_timings_seconds={};
events.evt_IDs={};

% and the markings (non spike)
infos.evt_timings_seconds={};
infos.evt_IDs={};


% start with empty offset, this will increase over time
offset=-1;
last_timestamp=0;

% now parse the events and adapt the offset if needed
for i=1:length(all_evts)
 if isempty(all_evts(i).value)
  % this is an empty marking, so check if 1 sec increase (with rounding tolerance) 
  if abs(all_evts(i).timestamp-last_timestamp-1)>0.0000000001
   % empty timestamp not in 1 sec interval, so assume new offset
   if offset<0
    % seems first item
    offset=all_evts(i).timestamp;
   else
    % adopt offset
    offset=all_evts(i).timestamp-last_timestamp+offset-1;
    warning('Break detected in the EDF(+) annotations. Attempting to change the offset to compensate.')
    fprintf('Note: This will not be 100%% exact (up to 1 second shift / break).\nIf precision is needed, export again as continous dataset (filling in gaps with 0).\nNew annotation offset: %0.4fsec\n',offset)
   end
   last_timestamp=all_evts(i).timestamp;
  else
   % is 1 second pulse, so keep as last legal
   last_timestamp=all_evts(i).timestamp;
  end
 else
  % not empty, so we treat as a label
  if ~isempty(regexpi(all_evts(i).value,'sp[ike]+'))
   % probably a spike marking
   events.evt_timings_seconds{end+1,1}=all_evts(i).timestamp-offset;
   events.evt_IDs{end+1,1}=all_evts(i).value;
  else
   % something else
   infos.evt_timings_seconds{end+1,1}=all_evts(i).timestamp-offset;
   infos.evt_IDs{end+1,1}=all_evts(i).value;
  end
 end
end

end

