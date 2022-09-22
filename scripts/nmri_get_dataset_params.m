function [ params ] = nmri_get_dataset_params( params, dataset_mapping )
% will modify the params according to dataset mapping
if ~isempty(dataset_mapping)
 if ~isstruct(dataset_mapping)
  % should be a file then
  if ~exist(dataset_mapping,'file')
   error(['could not find ' dataset_mapping])
  end
  % check if has params
  listOfVariables = who('-file', dataset_mapping);
  if ismember('params',listOfVariables)
   ds_params=load(dataset_mapping,'params');
   params=nf_update_struct(params,ds_params.params);
  end
 end
end

end

