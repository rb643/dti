function [Result] = runStats(Matrix, groupingVector, metric, covariate, nPerm, plotOn)
% function to run statistics on global graph metrics
% example: stats = runStats(nResult,p.Results.GroupMask,'dens',[],500,0); % to run stats for global network density
% run a standard anova for groups*metric and only runs post-hoc pairwise permutation is there is a main effect in the anova

%% get the number of unique group numbers
%% NB this assumes the groups are numbered sequentially
numgroups = length(unique(groupingVector));
[p,q] = meshgrid([1:numgroups], [1:numgroups]);
Result.pairs = [p(:) q(:)];

% split the groups, this might come in handy later
for i = 1:length(numgroups)
    mask = (groupingVector == i);
    group{i} = Matrix.(metric)(mask);
end

%%%%%% ANOVA %%%%%
data = Matrix.(metric);
[Result.p,Result.tbl,Result.stats] = anova1(data,groupingVector);

%%%%%% ANCOVA %%%%%
%[h,atab,ctab,stats] = aoctool(Matrix.(metric),covariate,groupingVector);


%%%%%% posthoc tests %%%%%
if Result.p < 0.05
    for i = 1:length(pairs)
        Result.permuted(i) = permutation_2tailed(group{Result.pairs(i,1)},group{Result.pairs(i,2)},nPerm);
        Result.meandiff(i) = mean(group{Result.pairs(i,1)}) - mean(group{Result.pairs(i,2)});
    end  
    Result.permuted';
    Result.meandiff';
end

Result;


end
