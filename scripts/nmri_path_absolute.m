function [outDat] = nmri_path_absolute(inDat)
%will take care of havin absolute paths from some relevant fields
%f=fieldnames(inDat);

% only check specific fields
f = cell(28, 1);
f{1,1} = 'raw_dataset';
f{2,1} = 'analysis_dir';
f{3,1} = 'dataset_mapping';
f{4,1} = 'layout';
f{5,1} = 'montages';
f{6,1} = 'elec_file';
f{7,1} = 'loaded_from';
f{8,1} = 'mri_dir';
f{9,1} = 'mri_ctf';
f{10,1} = 'mri_seg';
f{11,1} = 'suma_surface';
f{12,1} = 'clean_dataset';
f{13,1} = 'cleanICA_dataset';
f{14,1} = 'stats_dir';
f{15,1} = 'sensor_stats';
f{16,1} = 'stamp_dir';
f{17,1} = 'stamps';
f{18,1} = 'dws_filt_dataset';
f{19,1} = 'ICA_components';
f{20,1} = 'QCdir';
f{21,1} = 'suma_dir';
f{22,1} = 'fsdir';
f{23,1} = 'proc_dir';
f{24,1} = 'mri_T1';
f{25,1} = 'hdm_lead';
f{26,1} = 'electrodes_aligned';
f{27,1} = 'stats';
f{28,1} = 'SelectedTrials_file';

% make sure to only check existing fields
f=intersect(fieldnames(inDat),f);

outDat=inDat;

ft_dir=fileparts(which('ft_defaults'));

for i = 1:length(f)
 thisF=inDat.(f{i});
 if ischar(thisF) && ~strcmp(thisF(1),'/') && length(thisF)>4
  % we may want to root this field
  if exist(fullfile(pwd,thisF),'file')
   outDat.(f{i})=fullfile(pwd,thisF);
   % also check in Fieldtrip
  elseif exist(fullfile(ft_dir,'template','electrode',thisF),'file')
   outDat.(f{i})=fullfile(ft_dir,'template','electrode',thisF);
  end
 end
end


end

