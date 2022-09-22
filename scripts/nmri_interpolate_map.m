function [outmap, cfg] = nmri_interpolate_map(cfg,inmap)
%[outmap, cfg] = nmri_interpolate_map(cfg,inmap)
% 
% cfg              = config struct
%  .layout         = fieldtrip layout to use a reference
%  .neigh          = fieldtrip neighbourhood structure (will auto-determine
%                    if not set)
%  .dist           = Euclidian distances between neighbors (will auto-determine) 
%                   
% will assume that inmap is 1 value / layout channel/position and that the
% ordering of channels is concordant 
% for imporving processing speed neigh and dist may be re-used for the next
% run

if ~exist('cfg','var')
 error('Need a cfg variable')
end

if ~isfield(cfg,'layout')
 error('Need an layout info')
end

if ischar(cfg.layout)
 % try to load
 cfg = [];
 cfg.skipscale='yes';
 cfg.skipcomnt='yes';
 cfg.layout=ft_prepare_layout(cfg);
end

if ~isstruct(cfg.layout) || ~isfield(cfg.layout,'label') || ~isfield(cfg.layout,'pos')
 error('Layout is no Fieldtrip layout struct, cannot work with that')
end

if ~isfield(cfg,'neigh') || ~isstruct(cfg.neigh) 
 fprintf('Triangulating channel neighbours\n')
 cfg.neigh=ft_prepare_neighbours(struct('method','triangulation','layout',cfg.layout));
end

% make sure the labels match
if length(cfg.layout.label)~=length(cfg.neigh) || ~all(strcmp(cfg.layout.label,{cfg.neigh(:).label}'))
 error('Mismatch of layout labels and neighbours')
end

% calculate Euclidian distances (if not provided)
if ~isfield(cfg,'dist') || ~iscell(cfg.dist) || ~isfield(cfg,'neigh_idx') || ~iscell(cfg.neigh_idx) 
 fprintf('Calculating channel neighbours distances\n')
 cfg.dist=cell([length(cfg.layout.label),1]);
 cfg.neigh_idx=cell([length(cfg.layout.label),1]);
 for i=1:length(cfg.layout.label)
  cfg.neigh_idx{i}=cellfun(@(x) find(strcmp(x,cfg.layout.label)),cfg.neigh(i).neighblabel);
  cfg.dist{i}=pdist2(cfg.layout.pos(i,:),cfg.layout.pos(cfg.neigh_idx{i},:));
 end
end

if length(cfg.neigh_idx)~=length(cfg.dist)
 error('Mismatch of layout neigh distances and indices (should not happen)')
end


%% now find all the NaNs

% have maximum of N depth
loopC=1;

outmap=inmap;

% make sure our inmap is correctly ordered
if size(outmap,1)<size(outmap,2)
 outmap=outmap';
end

interpFields=isnan(outmap);
if size(outmap,1)>1 && size(outmap,2)>1
 % zero diagonal in 2D
 interpFields(eye(size(interpFields))==1)=0;
end
lastNaN=Inf;

while sum(interpFields(:))>0 && loopC<1000 && lastNaN>sum(interpFields(:))

 fprintf('Interpolation run N=%d, trying to replace %d NaNs\n',loopC,sum(interpFields(:)))
 lastNaN=sum(interpFields(:));
 thismap=outmap;
 % loop over max. 2 dimensions
 for x=1:size(outmap,1)
  for y=1:size(outmap,2)
   if interpFields(x,y) && (size(outmap,2)==1 || x~=y)
    if sum(~isnan(thismap(cfg.neigh_idx{x},y)))>0
     % now interpolate this value using the x / 1st dim as reference
     outmap(x,y)=nansum(thismap(cfg.neigh_idx{x},y).*(1./cfg.dist{x})')/sum(~isnan(thismap(cfg.neigh_idx{x},y)).*(1./cfg.dist{x})');
     % now check the inverse
    elseif size(thismap,2)>1 && sum(~isnan(thismap(x,cfg.neigh_idx{y})))>0
     % now interpolate this value using the y / 2nd dim as reference (for
     % 2D matricies)
     outmap(x,y)=nansum(thismap(x,cfg.neigh_idx{y}).*(1./cfg.dist{y})')/sum(~isnan(thismap(x,cfg.neigh_idx{y})).*(1./cfg.dist{y})'); 
    else
    end
   end
  end
 end
 loopC=loopC+1;
 interpFields=isnan(outmap);
 if size(outmap,1)>1 && size(outmap,2)>1
  % zero diagonal in 2D
  interpFields(eye(size(interpFields))==1)=0;
 end

end



end

