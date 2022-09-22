function [idx,p_age,p_sex] = nf_reg_match_ratio(regA,regB,ratio,strategy,p_cutoff,depth)
%[idx,p_age,p_sex] = nf_reg_match_ratio(regA,regB,ratio,strategy,p_cutoff,depth)
%   Function to automatically match two groups by age+sex 
%
% regA = regressor of reference group, i.e. the group that is matched to
% regB = regressor of variable group, i.e. the group that is matched /
%        reduced
% ratio = matching ration, 1:ratio (1 to ...), default 1
% strategy = 'sex>age', default
%               i.e. sex is matched first, then age
%            'age+sex', optional
%               sex and age are equally important, use an iteration
%               appraoch. I.e. we do a lot of random iterations/samplings
%               and pick the one that has the best match
% p_cutoff = p-value to accept for equality, default 0.5
% depth = vor internal usage (in the iter mode), do not change

if ~exist('depth','var')
 depth=1;
end

if ~exist('p_cutoff','var')
 p_cutoff=0.5;
end

if ~exist('strategy','var')
 strategy='sex>age';
end

if ~exist('ratio','var')
 ratio=1;
end

if strcmpi(strategy,'age+sex')
 % this is a special case, run iterations here
 fprintf('Doing 1000 random iterations to find the best variant...\n')
 aidx={};
 ap_age=[];
 ap_sex=[];
 % do 1000 iterations, since the 'age+sex' is rather random
 for ii=1:1000
  if mod(ii,100)==0
   fprintf('Iteration: %d\n',ii)
  end
  [aidx{ii},ap_age(ii),ap_sex(ii)] = nf_reg_match_ratio(regA,regB,ratio,'age+sex-iter',p_cutoff,depth);
 end
 % now find the best option
 fprintf('Mean Age p: %0.2f\n',mean(ap_age))
 fprintf('Mean Sex p: %0.2f\n\n',mean(ap_sex))
 
 [mp midx]=max((ap_age+ap_sex)/2);
 idx=aidx{midx};
 
else

 regBm={};

 s = RandStream('mt19937ar','Seed', seconds(round(milliseconds(second(datetime))*1000000)));

 % randomly order both
 available=1:length(regB.fname);
 available=available(randperm(s,length(available)));

 % randomly order both
 ref=1:length(regA.fname);
 ref=ref(randperm(s,length(ref)));


 sum_offset_age=0; % overall age offset
 sum_offset_female=0; % overall female excess


 if strcmpi(strategy,'sex>age') % 1st sex, then age

  for i=1:length(ref)
   % now search for matches
   cand=available(find(regB.grp_sex(available)==regA.grp_sex(ref(i))));

   if length(cand)<ratio
    fprintf('Found no or not enough match for case=%s and sex=%d, take any\n',regA.fname{ref(i)},regA.grp_sex(ref(i)))
    cand=available;
   end

   if length(cand)<ratio
    error('No candidates left, likely not enough cases in group B')
   end

   % now check the age distance
   found=0;
   while found<ratio
    d=regB.age(cand)-regA.age(ref(i))+sum_offset_age;
    [~,idx]=sort(abs(d));
    % take best match
    regBm=[regBm regB.fname(cand(idx(1)))];
    % remember our offset for next round
    sum_offset_age=sum_offset_age+d(idx(1));
    % and remove from candidates and available
    available(available==cand(idx(1)))=[];
    cand(idx(1))=[];
    found=found+1;
   end
  end

 elseif strcmpi(strategy,'age+sex-iter') %  age with priority, sex as minor criterion
  for i=1:length(ref)
   cand=available; % use all intially

   if length(cand)<ratio
    error('No candidates left, likely not enough cases in group B')
   end

   % now check the age distance
   found=0;
   while found<ratio
    d=regB.age(cand)-regA.age(ref(i))+sum_offset_age;
    % factor in female excess/miss
    [~,idx]=sort(abs(d)+(abs((regB.grp_sex(cand)-regA.grp_sex(ref(i))*-sum_offset_female))));
    % take best match
    regBm=[regBm regB.fname(cand(idx(1)))];
    % remember our offset for next round
    sum_offset_age=sum_offset_age+d(idx(1));
    sum_offset_female=regA.grp_sex(ref(i))-regB.grp_sex(cand(idx(1)));
    % and remove from candidates and available
    available(available==cand(idx(1)))=[];
    cand(idx(1))=[];
    found=found+1;
   end
  end
 else
  error(['Unknown mataching strategy ' strategy])
 end


 % return the original indices, sorted
 idx=sort(cellfun(@(x) find(strcmp(regB.fname,x)),regBm));

end


% now check if the distribution is okay
if ~strcmpi(strategy,'age+sex-iter')
 fprintf('Age:\n%0.1f +/- %0.1f for reference\n%0.1f +/- %0.1f for matched group\n',mean(regA.age),std(regA.age),mean(regB.age(idx)),std(regB.age(idx)))
 fprintf('Min: %d, Max: %d for reference\nMin: %d, Max: %d for matched group\n\n',min(regA.age),max(regA.age),min(regB.age(idx)),max(regB.age(idx)))
end

[h,p_age,ks2stat] = kstest2(regA.age,regB.age(idx));
%[h,p_age,tstat] = ttest2(regA.age,regB.age(idx));


if ~strcmpi(strategy,'age+sex-iter')
 fprintf('Kolmogorov-Smirnov test for coming from the same distribution is:\n');
 %fprintf('T-test for coming from the same distribution is:\n');

 fprintf('p: %0.4f\n\n',p_age);

 fprintf('Sex:\n%d females, %d males for reference\n%d females, %d males for matched group\n',sum(regA.grp_sex==1),sum(regA.grp_sex==0),sum(regB.grp_sex(idx)==1),sum(regB.grp_sex(idx)==0))
 fprintf('%0.1f%% females, %0.1f%% males for reference\n%0.1f%% females, %0.1f%% males for matched group\n',sum(regA.grp_sex==1)*100/length(regA.grp_sex),sum(regA.grp_sex==0)*100/length(regA.grp_sex),...
 sum(regB.grp_sex(idx)==1)*100/length(regB.grp_sex(idx)),sum(regB.grp_sex(idx)==0)*100/length(regB.grp_sex(idx)))
end

%[h,p_sex,ks2stat] = kstest2(regA.grp_sex,regB.grp_sex(idx));
[tbl,ch2,p_sex] = crosstab([regA.grp_sex ; regB.grp_sex(idx)], [ones(length(regA.grp_sex),1); zeros(length(idx),1)]);

if ~strcmpi(strategy,'age+sex-iter')
 fprintf('Chi-Square test for coming from the same distribution is:\n');
 fprintf('p: %0.4f\n\n',p_sex);
end


p=min(p_sex,p_age);

if p>p_cutoff
 if ~strcmpi(strategy,'age+sex-iter')
  fprintf('This seems okay!\n\n')
 end
else
 depth=depth+1;
 if depth>10
  if ~strcmpi(strategy,'age+sex-iter')
   fprintf('Even after 10 runs, we did not succseed here. Abortion.\n\n')
  end
 else
  if ~strcmpi(strategy,'age+sex-iter')
   fprintf('This seems NOT okay! Calling the fitting again (new depth=%d)\n\n',depth)
  end
  % recurse
  [idx,p_age,p_sex] = nf_reg_match_ratio(regA,regB,ratio,strategy,p_cutoff,depth);
 end
end


end

