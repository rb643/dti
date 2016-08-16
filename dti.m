function [Result wResult pResult faResult] = dti(Matrix,varargin)
% 
% run as: [nResult wResult pResult faResult] = dti('Matrix',adjacencymatrix);
% 
% ---------------------------- INPUT ----------------------------
%
% Matrix            -       4D connectivity matrix(node*node*measure*subject)
%
% --------------------- OPTIONAL ARGUMENTS ----------------------
%
% regionDescriptions-       A cell structure containing names of all nodes in your matrix (at this point its only used for plotting)
% nRand             -       number of randomization used in random networks for normalization
% PlotAverage       -       Logical to set if we want to plot average
% PlotAll           -       Logical to set if we want to plot the output of some graph metrics
% prev              -       prevalence for prevalence weighted matrices
% metric            -       metric that we want to use to construct matrices (defaults to number of streamlines)
% groups            -       1D vector with grouplabels
%
% ---------------------------- OUTPUT ----------------------------
% Result            -       Graph metrics from binary matrices including all NOS's >0
% wResult           -       Graph metrics from weighted matrices
% pResult           -       Graph metrics from matrices weighted based on group prevalence
% TODO:
% nResult           -       Graph metrics from matrices corrected with ROI*ROI volume

% TODO:
% include subject mask to tell which subject to include (although this is probably better to do beforehand...)
% include ROI volume vector for volume weighting/correction


%% check all the inputs and if they do not exist then revert to default settings
% set the larger default in case they are not specified
groupdefault = [1 2 3 0	2 0	0	0	3	0	2	3	0	1	1	0	0	0	1	0	2	2	2	1	2	3	3	3	3	2	1	3	2	3	0	0	2	0	1	2	0	0	0	1	1	1	3	2	3	2	1	3	1	1	1	2	1	1	2	1	1	2	1	1	1	2	1	2	1	0	2	1	3	1	1	0	0	3	0	0	0	0	2	1	2	3	2	2	1	1	6	1	3	3	2	3	3]';
regionsdefault = num2str([1:82]);
% input parsing settings
p = inputParser;
p.CaseSensitive = true;
p.Parameters;
p.Results;
p.KeepUnmatched = true;
% set the desired and optional input arguments
addRequired(p,'Matrix',@isnumeric);
addOptional(p,'regionDescriptions',regionsdefault,@iscell);
addOptional(p,'nRand',10,@isnumeric);
addOptional(p,'PlotAverage',0,@isnumeric);
addOptional(p,'PlotAll',0,@isnumeric);
addOptional(p,'prev',0.75,@isnumeric);
addOptional(p,'metric',1,@isnumeric);
addOptional(p,'groups',groupdefault, @isnumeric);
% parse the input
parse(p,varargin{:});
% then set get all the inputs out of a structure
nRand = p.Results.nRand; metric = p.Results.metric; regionDescriptions = p.Results.regionDescriptions; PlotAverage = p.Results.PlotAverage; PlotAll = p.Results.PlotAll; prev = p.Results.prev; groups = p.Results.groups; connectivity = p.Results.Matrix;

%% And here we go with some actual analysis
nSubjects = size(connectivity,4);
% Extract NOS matrix for every subject and remove singleton dimensions
NOS_w = squeeze(connectivity(:,:,metric,:));
% Binarise NOS matrix for every subject
NOS_b = double(NOS_w > 0);

%% Inserting group membership and demographic and cognitive information
Group0 = find(groups == 0); Group1 = find(groups == 1); Group2 = find(groups == 2); Group3 = find(groups == 3); Group4 = find(groups == 4);

%% Quality control
% Compute prevalence matrix of entire sample which will be used later.
P = mean(connectivity(:,:,metric,:) > 0, 4);
% Create a vector of average NOS over all connections of a subject
for i = 1:nSubjects
    NOS = connectivity(:,:,metric,i);
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
    AvgPrevMissingCon(i,1) = mean(P(C == 0));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Unweighted NOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate degree, density, length, and clustering for each subject
% Everything here is done on a single-level, only use prevalence-thresholds
% when looking at group average connectomes
for i = 1:nSubjects
    A = double(connectivity(:,:,metric,i) > 0);
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
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Binary Networks based on NOS','HorizontalAlignment','center','VerticalAlignment', 'top');
    
    figure;
    subplot(4,1,1); bar(mean(Result.deg(groups == 0,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionDescriptions)),'XLim',[0 (length(regionDescriptions)+1)],'XTickLabel',regionDescriptions, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 0'); title('Degree');
    subplot(4,1,2); bar(mean(Result.deg(groups == 1,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionDescriptions)),'XLim',[0 (length(regionDescriptions)+1)],'XTickLabel',regionDescriptions, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 1');
    subplot(4,1,3); bar(mean(Result.deg(groups == 2,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionDescriptions)),'XLim',[0 (length(regionDescriptions)+1)],'XTickLabel',regionDescriptions, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 2');
    subplot(4,1,4); bar(mean(Result.deg(groups == 3,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionDescriptions)),'XLim',[0 (length(regionDescriptions)+1)],'XTickLabel',regionDescriptions, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Weighted NOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Weighted networks
for i = 1:nSubjects
    A = squeeze(connectivity(:,:,metric,i));
    % create non-normalized output
    wResult.deg(i,:) = degrees_und(A); % degree
    wResult.dens(i,:) = density_und(A); %density
    wResult.cpl(i,:) = charpath(distance_wei(A)); %characteristic path length
    wResult.clust(i,:) = mean(clustering_coef_wu(A)); %clustering coefficient
    wResult.trans(i,:) = transitivity_wu(A); %transitivity
end

if PlotAll == 1
    figure;
    subplot(2,2,1); boxplot(wResult.dens,groups,'colorgroup',groups);title('density');
    subplot(2,2,2); boxplot(wResult.cpl,groups,'colorgroup',groups);title('characteristic path length');
    subplot(2,2,3); boxplot(wResult.clust,groups,'colorgroup',groups);title('clustering');
    subplot(2,2,4); boxplot(wResult.trans,groups,'colorgroup',groups);title('transitivity');
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Weighted Networks','HorizontalAlignment','center','VerticalAlignment', 'top');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Unweighted by prevalence %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Network with NOS in >40% of subject
PrevalenceMask = squeeze(double(connectivity(:,:,1,:)>0));
PrevalenceMask = mean(PrevalenceMask,3);
PrevalenceMask = double(PrevalenceMask>prev);

for i = 1:nSubjects
    A = double(connectivity(:,:,metric,i) > 0);
    A = A.*PrevalenceMask;
    
    % create non-normalized output
    pResult.deg(i,:) = degrees_und(A); % degree
    pResult.dens(i,:) = density_und(A); %density
    pResult.cpl(i,:) = charpath(distance_bin(A)); %characteristic path length
    pResult.clust(i,:) = mean(clustering_coef_bu(A)); %clustering coefficient
    pResult.rich(i,:) = rich_club_bu(A,35); %rich-club coefficient
    pResult.trans(i,:) = transitivity_bu(A); %transitivity
    pResult.assor(i,:) = assortativity_bin(A,0); %assortativity
    pResult.effN(i,:) = efficiency_nodal(A); %efficiency
    
    % create a random matrix for normalization
    R = randmio_und(A,nRand);
    Crand = mean(clustering_coef_bu(R)); % get random clustering
    Lrand = charpath(distance_bin(R)); % get path length from random matrix
    pResult.Sigma(i,:)=(pResult.clust(i,:)./Crand)./(pResult.cpl(i,:)./Lrand); % get small world coefficient
    
    for iR = 1:nRand
        TempClub(iR,:,:) = randmio_und_connected(A, nRand);% create a random matrix from original
        RichRand(iR,:) = rich_club_bu(squeeze(TempClub(iR,:,:)),35);
    end
    RichRand = mean(RichRand,1);
    
    pResult.Norm.rich(i,:) = pResult.rich(i,:)./RichRand;
end

if PlotAll == 1
    figure;
    subplot(2,2,1); boxplot(pResult.dens,groups,'colorgroup',groups);title('density');
    subplot(2,2,2); boxplot(pResult.cpl,groups,'colorgroup',groups);title('characteristic path length');
    subplot(2,2,3); boxplot(pResult.clust,groups,'colorgroup',groups);title('clustering');
    subplot(2,2,4); boxplot(pResult.Sigma,groups,'colorgroup',groups);title('small-world coefficient');
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Weighted according to Prevalence','HorizontalAlignment','center','VerticalAlignment', 'top');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Weighted by NOS corrected for VOI %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Weighted by FA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:nSubjects
    A = squeeze(connectivity(:,:,3,i));
    % create non-normalized output
    faResult.deg(i,:) = degrees_und(A); % degree
    faResult.dens(i,:) = density_und(A); %density
    faResult.cpl(i,:) = charpath(distance_wei(A)); %characteristic path length
    faResult.clust(i,:) = mean(clustering_coef_wu(A)); %clustering coefficient
    faResult.trans(i,:) = transitivity_wu(A); %transitivity
end

if PlotAll == 1
    figure;
    subplot(2,2,1); boxplot(faResult.dens,groups,'colorgroup',groups);title('density');
    subplot(2,2,2); boxplot(faResult.cpl,groups,'colorgroup',groups);title('characteristic path length');
    subplot(2,2,3); boxplot(faResult.clust,groups,'colorgroup',groups);title('clustering');
    subplot(2,2,4); boxplot(faResult.trans,groups,'colorgroup',groups);title('transitivity');
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Weighted Networks by FA','HorizontalAlignment','center','VerticalAlignment', 'top');
end

end
