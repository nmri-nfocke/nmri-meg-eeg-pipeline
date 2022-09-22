function nmri_combine_exports(inputExports,pathOut)
%nmri_combine_exports(inputExports,pathOut)
%   Function to combine multiple (>=2) export EEG/MEG export directories
%   into a single one
%
% inputExports = cell array of paths to exports to combine
% pathOut      = new output path to be generated


if ~exist('inputExports','var') || ~iscell(inputExports)
 error('Need a cell array of exports paths')
end


if ~exist('pathOut','var') || ~ischar(pathOut)
 error('Need a char string path for unified export ')
end

if length(inputExports)<2
 error('Joining needs at least 2 export paths')
end

% now list the 1st path
files=dir(inputExports{1});

% we need to make sure to have a .mat file 1st to allow full parsing


% and loop these
newN=-1;
for i=1:length(files)
 % skip .
 if strcmp(files(i).name(1),'.')
  continue
 end
 
 [~,fi,ext]=fileparts(files(i).name);
 
 % now check for the same file in other paths
 matchPattern=regexprep(files(i).name,'_N[0-9]*','_N[0-9]*')
 foundMatch={fullfile(files(i).folder,files(i).name)};
 for ii=2:length(inputExports)
  filesM=dir(inputExports{ii});
  fidx=find(~cellfun(@isempty,regexp({filesM.name},matchPattern)));
  if length(fidx)==1
   % single match found
   foundMatch=[foundMatch;{fullfile(filesM(fidx).folder,filesM(fidx).name)}];
  else
   fprintf('Could not find a match for %s in %s\n',matchPattern,inputExports{ii})
  end
 end
 
 if length(foundMatch)~=length(inputExports)
  warning(sprintf('Incomplete file sets for %s. Joined export may be defect.',matchPattern))
  continue
 end
 
 % if we are here, we have a set of files, so try to make some sense of it
 if strcmp(ext,'.mat')
  % seems a .mat file, attempt to combine 
  combMat=[];
  for ii=1:length(foundMatch)
   thisMat=load(foundMatch{ii});
   fields=fieldnames(thisMat);
   for fi=1:length(fields)
    % check active_hdm_class
    if strcmp(fields{fi},'active_hdm_class')
     if isfield(combMat,'active_hdm_class')
      if ~strcmp(combMat.active_hdm_class,thisMat.active_hdm_class)
       error('Mismatch of active_hdm_class detected. Fatal.')
      end
     else
      combMat.active_hdm_class=thisMat.active_hdm_class;
     end
    elseif strcmp(fields{fi},'all_msk')
     % the logical mask
     if isfield(combMat,'all_msk')
      % and combine
      combMat.all_msk=combMat.all_msk&thisMat.all_msk;
     else
      combMat.all_msk=thisMat.all_msk;
     end
     
    else
     % all other fields should be concatable
     if isfield(combMat,fields{fi})
      % vert cat
      combMat.(fields{fi})=[combMat.(fields{fi});thisMat.(fields{fi})];
     else
      % make
      combMat.(fields{fi})=thisMat.(fields{fi});
     end
    end
   end
  end
  
  % determine N
  fields=fieldnames(combMat);

  for fi=1:length(fields)
   if ~strcmp(fields{fi},'active_hdm_class') && ~strcmp(fields{fi},'all_msk')
    if newN<0
     % set newN
     newN=length(combMat.(fields{fi}));
    else
     if length(combMat.(fields{fi}))~=newN
      error('Mismatch of N detected, should not happen. Fatal.')
     end
    end
   end
   
   
   % now write out
   outfile=fullfile(pathOut,regexprep(files(i).name,'_N[0-9]*',['_N' num2str(newN)]));
   if ~exist(pathOut,'dir')
    mkdir(pathOut)
   end
   save(outfile,'-struct','combMat') 
  end
  
 elseif strcmp(ext,'.csv') || strcmp(files(i).name(end-3:end),'_log')
  % write sequentially
  combTxt='';
  for ii=1:length(foundMatch)
   thisTxt=fileread(foundMatch{ii});
   combTxt=[combTxt,thisTxt];
  end
  
  % determine N
  thisN=length(regexp(combTxt,'[\n]'));
  if newN<0
   % set newN
   newN=thisN;
  else
   if thisN~=newN
    %error('Mismatch of N detected, should not happen. Fatal.')
   end
  end
  
  % now write out
  outfile=fullfile(pathOut,regexprep(files(i).name,'_N[0-9]*',['_N' num2str(newN)]));
  fid=fopen(outfile,'w');
  fprintf(fid,'%s',combTxt);
  fclose(fid);
  
  
 end
end
