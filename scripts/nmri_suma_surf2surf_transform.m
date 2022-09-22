function [surfB_values, surfB_dist, nn_vertices, nn_weights ]=nmri_suma_surf2surf_transform(surfA, surfB, surfA_values, nn_vertices, nn_weights )
%nmri_suma_surf2surf_transform(surfA, surfB, surfA_values, nn_vertices, nn_weights )
%
% Transforms a vertex-based parameter/value to a new surface based in the
% shortest distances between vertices (A and B). 
% Note: Does not do any registration between surfaces, need to be done
% beforehands.
% also note that SUMA (ld.??) surface coord. are NOT in register with
% Freesurfer sphere and need to be mapped differently
%
% nn_weights /_vertices can be re-used for speed if using the same mapping
% again


% do some safety checks
if (size(surfA.pos,1)~=size(surfA_values,1))
 error('Mismatch of surfA vertices and surfA values')
end

if  ~exist('nn_vertices','var') && nargout>2
 % save the weigths
 nn_vertices=cell(length(surfB.pos),1);
 nn_weights=cell(length(surfB.pos),1);
end

% low loop over each target vertex and create mapping
surfB_values=zeros(length(surfB.pos),1);
surfB_dist=zeros(length(surfB.pos),1);
for i=1:length(surfB.pos)
 if ~exist('nn_vertices','var') || isempty(nn_vertices{i})
  % not provided - calculate
  xyz=surfB.pos(i,1:3);
  % restrict to same anot-key and hemisphere
  this_annot=find((surfA.annot==surfB.annot(i,1))&(surfA.hemi==surfB.hemi(i,1)));
  d=pdist2(surfA.pos(this_annot,1:3),xyz);
  [ds,iA]=sort(d); %sort and keep indicies
  surfB_dist(i,1)=ds(1);
  % now find the original vertex
  ref_vertex=this_annot(iA(1));
  % now find legitimate faces of closest match
  [row, ~]=find(surfA.tri==ref_vertex); 
  n_vert=unique(surfA.tri(row)); % these are connected neighbours
  % use only same annot
  n_vert=intersect(n_vert,this_annot);
  % check if we have at least 3
  c=0;
  while length(n_vert)<6 & ~isempty(n_vert) & c<10
   % then check neighbour's neighbours also
   nn_vert=[];
   for ii=1:length(n_vert)
    [row, ~]=find(surfA.tri==n_vert(ii));
    nn_vert=unique([surfA.tri(row);nn_vert]);
   end
   % use only same annot
   nn_vert=intersect(nn_vert,this_annot);
   n_vert=unique([n_vert;nn_vert]);
   c=c+1;
  end
  
  % get distance again, using double for precision
  d=pdist2(double(surfA.pos(n_vert,1:3)),double(xyz));
  [~,n_i]=sort(d); %sort them by distance
  % and take closest 9, if we have them
  nn_i=[];
  for ii=1:9
   if length(n_i)>=ii
    nn_i=[nn_i n_i(ii)];
    %avoid divide by 0 and reduce perfect match impact
    if (d(n_i(ii))<0.1)
     d(n_i(ii))=0.1;
    end
   end
  end
  % now make weights
  n_weights=1./(d(nn_i).^1.5);
  w_sum=sum(reshape(n_weights,numel(n_weights),1));
  % standardize to 1
  n_weights=n_weights/w_sum;
  % get the values
  this_vals=surfA_values(n_vert(nn_i),1);
  % deal with NaN
  if sum(~isnan(this_vals))==0
   surfB_values(i,1)=NaN;
  else
   % not all are Nan, use them
   t_weights=n_weights(~isnan(this_vals));
   t_sum=sum(reshape(t_weights,numel(t_weights),1));
   % standardize to 1
   t_weights=t_weights/t_sum;
   % remove NaN
   this_vals(isnan(this_vals))=[];
   this_vals=this_vals.*t_weights;
   surfB_values(i,1)=sum(reshape(this_vals,numel(this_vals),1));
  end
  if nargout>2
   % remember if wanted
   nn_vertices{i}=n_vert(nn_i);
   nn_weights{i}=n_weights;
  end
 else
  % neighbours and weights provided
  % get the values
  this_vals=surfA_values(nn_vertices{i},1);
  % deal with NaN
  if sum(~isnan(this_vals))==0
   surfB_values(i,1)=NaN;
  else
   % not all are Nan, use them
   t_weights=nn_weights{i}(~isnan(this_vals));
   t_sum=sum(reshape(t_weights,numel(t_weights),1));
   % standardize to 1
   t_weights=t_weights/t_sum;
   % remove NaN
   this_vals(isnan(this_vals))=[];
   this_vals=this_vals.*t_weights;
   surfB_values(i,1)=sum(reshape(this_vals,numel(this_vals),1));
  end 
 end
end


end
