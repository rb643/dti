function [Result] = eli(Matrix, GroupMask, nRand, PlotAverage, PlotAll)
%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%
%% Matrix            -       4D connectivity matrix(node*node*measure*subject)
%% GroupMask         -       1D vector with grouplabels
%% UseMask           -       1D logical with which subjects to include, if not supplied defaults to all
%% nRand             -       number of randomization used in random networks for normalization
%% PlotAverage       -       Logical to set if we want to plot average
%% PlotAverage       -       Logical to set if we want to plot the output of some graph metrics 


if isempty(GroupMask)
    % Assign default group membership
    groups = [1 2 3 0	2 0	0	0	3	0	2	3	0	1	1	0	0	0	1	0	2	2	2	1	2	3	3	3	3	2	1	3	2	3	0	0	2	0	1	2	0	0	0	1	1	1	3	2	3	2	1	3	1	1	1	2	1	1	2	1	1	2	1	1	1	2	1	2	1	0	2	1	3	1	1	0	0	3	0	0	0	0	2	1	2	3	2	2	1	1	6	1	3	3	2	3	3]';
else
    groups = GroupMask;
end

if isempty(nRand); nRand = 10; end
if isempty(PlotAverage); PlotAverage = 1; end

connectivity = Matrix;
nSubjects = size(connectivity,4);
% Extract NOS matrix for every subject and remove singleton dimensions
NOS_w = squeeze(connectivity(:,:,1,:));

% Binarise NOS matrix for every subject
NOS_b = double(NOS_w > 0);

%% Inserting group membership and demographic and cognitive information


Group0 = find(groups == 0);
Group1 = find(groups == 1);
Group2 = find(groups == 2);
Group3 = find(groups == 3);
Group4 = find(groups == 4);

% Demographics and cognitive variables
Age = [75	78	74	69	81	78	78	75	67	71	74	87	79	81	81	71	73	73	78	76	79	64	75	71	69	69	73	81	79	70	66	78	84	79	81	82	81	79	81	83	72	63	62	74	70	78	79	69	77	62	87	80	78	78	77	76	65	84	73	88	79	80	81	80	77	87	71	89	78	71	80	76	70	62	73	80	78	70	75	66	80	75	79	89	81	88	88	67	66	74	73	77	76	80]';
Gender = [1	1	2	2	2	1	2	1	1	2	2	2	1	1	2	1	2	1	1	1	2	1	1	1	1	1	1	1	1	1	1	1	1	1	1	2	2	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	2	1	1	1	1	1	1	2	1	1	2	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	2	2	1	2	2	1	1	1	1	1	1	2	1	1]';
MMSE = [30	20	28	30	28	23	30	25	28	28	21	15	28	29	28	27	30	27	20	28	18	16	21	21	21	26	26	12	25	28	21	29	28	28	30	18	25	29	29	30	23	26	19	26	22	27	19	22	25	23	15	17	23	24	22	25	23	22	29	22	19	25	28	26	23	17	28	22	24	23	18	23	29	30	23	30	30	29	30	26	19	18	28	19	19	25	23	20	25	26	25	16	19	16]';

%% Quality control
% Compute prevalence matrix of entire sample which will be used later.
P = mean(connectivity(:,:,1,:) > 0, 4);

% Create a vector of average NOS over all connections of a subject
for i = 1:nSubjects
    NOS = connectivity(:,:,1,i);
    average_NOS(i,1) = mean(NOS(:)); %% should this be the mean of the entire ? Or would it be better to first tril the matrix?
    % Create a vector of average FA over all connections of a subject
    FA = connectivity(:,:,3,i);
    average_FA(i,1) = mean(FA(:));
    % Create average prevalence of connections of a subject
    W = connectivity(:,:,1,i);
    B = double(W > 0);
    AvgPrevCon(i,1) = mean(P(B ==1));
    % Create average prevalence of missing connections of a subject
    W = connectivity(:,:,1,i);
    B = double(W > 0);
    C = B + eye(size(B));
    AvgPrevMissingCon(i,1) = mean(P(C == 0)); %why need 1 here and not for AvgPrevCon?
end

% Store all QC variables into a QC matrix
matrix_QC = [average_NOS average_FA AvgPrevCon AvgPrevMissingCon];

% Plot output for QC
if PlotAverage == 1
    figure;
    subplot(2,2,1); boxplot(matrix_QC(:,1),groups,'colorgroup',groups);title('NOS');
    subplot(2,2,2); boxplot(matrix_QC(:,2),groups,'colorgroup',groups);title('FA');
    subplot(2,2,3); boxplot(matrix_QC(:,3),groups,'colorgroup',groups);title('Avg Prev of Connections');
    subplot(2,2,4); boxplot(matrix_QC(:,4),groups,'colorgroup',groups);title('Avg Prev of Missing Connections');
end

%% Global metrics
% Create prevalence matrix of the sample
P_Sample = mean(connectivity(:,:,1,:) > 0, 4);
A_Sample = double(P>=0.6);

P_0= mean(connectivity(:,:,1,Group0) > 0, 4);
A_0 = double(P_0>=0.5);

P_1 = mean(connectivity(:,:,1,Group1) > 0, 4);
A_1 = double(P_1>=0.5);

P_2= mean(connectivity(:,:,1,Group2) > 0, 4);
A_2 = double(P_2>=0.5);

P_3 = mean(connectivity(:,:,1,Group3) > 0, 4);
A_3= double(P_3>=0.5);

if PlotAverage == 1
    figure;
    subplot(2,2,1);
    imagesc(A_0); title('0');
    subplot(2,2,2);
    imagesc(A_1); title('1');
    subplot(2,2,3);
    imagesc(A_2); title('2');
    subplot(2,2,4);
    imagesc(A_3); title('3');
end

% Create degree distribution
degrees_0 = sum(A_0,2); degrees_1 = sum(A_1,2); degrees_2 = sum(A_2,2); degrees_3 = sum(A_3,2);

if PlotAverage == 1
    figure;
    subplot(2,2,1);
    hist(degrees_0); title('0');  xlabel('deg of prevalence > 0.4')
    subplot(2,2,2);
    hist(degrees_1); title('1');  xlabel('deg of prevalence > 0.4')
    subplot(2,2,3);
    hist(degrees_2); title('2');  xlabel('deg of prevalence > 0.4')
    subplot(2,2,4);
    hist(degrees_3); title('3');  xlabel('deg of prevalence > 0.4')
end

%% Calculate degree, density, length, and clustering for each subject
% Everything here is done on a single-level, only use prevalence-thresholds
% when looking at group average connectomes
for i = 1:nSubjects
    A = double(connectivity(:,:,1,i) > 0);
    % create non-normalized output
    Result.deg(i,:) = degrees_und(A); % degree
    Result.dens(i,:) = density_und(A); %density
    Result.cpl(i,:) = charpath(distance_bin(A)); %characteristic path length
    Result.clust(i,:) = mean(clustering_coef_bu(A)); %clustering coefficient
    Result.rich(i,:) = rich_club_bu(A,35); %rich-club coefficient
    Result.trans(i,:) = transitivity_bu(A); %transitivity
    Result.assor(i,:) = assortativity_bin(A,0); %assortativity
    Result.effN(i,:) = efficiency_nodal(A); %efficiency
    
     % create a random matrix for normalization
    R = randmio_und(A,nRand);
        Crand = mean(clustering_coef_bu(R)); % get random clustering
        Lrand = charpath(distance_bin(R)); % get path length from random matrix
        Result.Sigma(i,:)=(Result.clust(i,:)./Crand)./(Result.cpl(i,:)./Lrand); % get small world coefficient
    cpl_random(i,:) = mean(squareform(distance_bin(R)));
    clust_random(i,:) = mean(clustering_coef_bu(R));
    for iR = 1:nRand
        TempClub(iR,:,:) = randmio_und_connected(A, nRand);% create a random matrix from original
        RichRand(iR,:) = rich_club_bu(squeeze(TempClub(iR,:,:)),35);
    end
    RichRand = mean(RichRand,1);
    
    Result.Norm.rich(i,:) = Result.rich(i,:)./RichRand;
end

if PlotAll == 1
    figure;
    subplot(2,2,1); boxplot(Result.dens,groups,'colorgroup',groups);title('density');
    subplot(2,2,2); boxplot(Result.cpl,groups,'colorgroup',groups);title('characteristic path length');
    subplot(2,2,3); boxplot(Result.clust,groups,'colorgroup',groups);title('clustering');
    subplot(2,2,4); boxplot(Result.Sigma,groups,'colorgroup',groups);title('small-world coefficient');
end

% % generating the mean normalised CPL and CLUS for every subject
% cpl_random = mean(cpl_random,2); clust_random = mean(clust_random,2);
% cpl_norm = cpl./cpl_random; clust_norm = clust./clust_random;
% 
% % Plot boxplot of normalised measures of path length and clustering across the groups
% figure(8)
% subplot(2,2,1);
% boxplot(cpl_norm, groups)
% title('Norm Path Length')
% subplot(2,2,2);
% boxplot(clust_norm, groups)
% title('Norm Clustering')
% subplot(2,2,3);
% boxplot(cpl, groups)
% title('Path Length')
% subplot(2,2,4);
% boxplot(clust, groups)
% title('Clustering')
% 
% %% Modules and hubs
% % Richard: I think we can use the hub-section of your recent manuscript for
% % this too. Not very confident about what I have written for hubs here so
% % feel free to revamp it :) This was just based off some of the examples in
% % the MVH's summer school book.
% 
% % Divide the nodes of the group-averaged connectome maps (0.4%) into
% % modules
% [M_0 Q_0] = modularity_und(A_0);
% [M_1 Q_1] = modularity_und(A_1);
% [M_2 Q_2] = modularity_und(A_2);
% [M_3 Q_3] = modularity_und(A_3);
% 
% %% Visualise modules across groups
% [~,I_0] = sort(M_0);
% [~,I_1] = sort(M_1);
% [~,I_2] = sort(M_2);
% [~,I_3] = sort(M_3);
% 
% figure(10)
% subplot(2,2,1);
% imagesc(A_0(I_0,I_0)); colorbar, title('0')
% subplot(2,2,2);
% imagesc(A_1(I_1,I_1)); colorbar;title('1')
% subplot(2,2,3);
% imagesc(A_2(I_2,I_2)); colorbar;title('2')
% subplot(2,2,4);
% imagesc(A_3(I_3,I_3)); colorbar;title('3')
% 
% 
% % How can we compare hubs across groups? I think the reviewer raised a
% % simiar point in your manuscript! Here I am trying to extract Q.
% % Q is the modularity of the network and indicates to what extent the
% % nework can be subdivided into separate modules.
% 
% % Incompleted.
% 
% 
% %% Hubs
% % Find the 12 highest degree of the binary group averaged connectome map A.
% % I guess this is not very relevant to the group pcomparisons as it across
% % the entire sample
% 
% degree_of_node_across_sample = sum(A,2);
% [~,I_degree_across_sample] = sort(degree_of_node_across_sample,'descend');
% regionDescriptions(I_degree_across_sample(1:12,:))
% 
% 
% % How to generate group averaged hub maps?
% degree_of_node_across_1 = sum(A_1,2);
% [~,I_degree_across_sample_1] = sort(degree_of_node_across_1,'descend');
% regionDescriptions(I_degree_across_sample_1(1:5))
% 
% degree_of_node_across_0 = sum(A_0,2);
% [~,I_degree_across_sample_0] = sort(degree_of_node_across_0,'descend');
% regionDescriptions(I_degree_across_sample_0(1:5))
% 
% degree_of_node_across_2 = sum(A_2,2);
% [~,I_degree_across_sample_2] = sort(degree_of_node_across_2,'descend');
% regionDescriptions(I_degree_across_sample_2(1:5))
% 
% degree_of_node_across_3 = sum(A_3,2);
% [~,I_degree_across_sample_3] = sort(degree_of_node_across_3,'descend');
% regionDescriptions(I_degree_across_sample_3(1:5))
% 
% % Now that the hubs are identified, what are the hub differences between
% % the groups? Compare strength of hubs in A vs B.
% % Or, take the hubs in HC, extrapolate hubs to other groups, compare metrics
% 
% % Plot them in BrainNetViewer using MNI coordinates
% 
% % Compute betweenness centrality of every node
% betweenness_bin(A);
% 
% %%  Rich Club
% 
% % As above, this is what I have written. Very elementary stuff for the
% % analyses of rich club networks. Feel free to expand.
% 
% R_0= rich_club_bu(A_0);
% R_1 = rich_club_bu(A_1);
% figure; plot(1:numel(R), R, '-o')
% 
% % Generate 100 randomised networks and save their rich club ceofficient
% % vectors as rows of matrix Rrandom. So it should have 100 rows for each of
% % the network and numel(R) columns.
% 
% Rrandom_0 = [zeros(numel(100),numel(R_0))];
% %The R here is the same as the R defined previously in the whole sample.
% 
% for i = 1:100
%     B_0 = randmio_und(A_0,5);
%     Rrandom_0(i,:) = rich_club_bu(B_0);
% end
% 
% Rrandom_1 = [zeros(numel(100),numel(R))];
% %The R here is the same as the R defined previously in the whole sample.
% 
% for i = 1:100
%     B_1 = randmio_und(A_1,5);
%     Rrandom_1(i,:) = rich_club_bu(B_1);
% end
% 
% % Plot rich club and corresponding random curves
% figure; plot(1:numel(R_HC), R_HC, '-o', 1:numel(R_HC), mean(Rrandom_HC, 1), '-o',1:numel(R_AD), R_AD, '-o', 1:numel(R_AD), mean(Rrandom_AD, 1), '-o')
% figure; plot(1:numel(R_HC), R_HC, '-o', 1:numel(R_AD), R_AD, '-o')
% 
% % Form vector hubs revealed earlier
% % Compute density
% % Create logical matrix of rich club, feeder, and local
% % Do weighted rich club and statistical comparisons between groups
% 
% %% Weighted graph metrics
% % Richard: All the analyses have been done on binarised networks so far, I'd like to
% % perform analyses on the weighted matrices too!

end

