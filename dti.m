function [nResult wResult pResult p] = dti(Matrix,varargin)
%
% For example:
% [nResult wResult pResult p] = dti('Matrix',adjacencymatrix);
% [nResult wResult pResult p] = dti('Matrix',adjacencymatrix,'regionLabels',regions,'prev',0.75,'PlotLocal',1,'PlotGlobal',...
%                                   1,'nos',10, 'subjectmask',include, 'groups', groups);
%
% ---------------------------- INPUT ----------------------------
%
% Matrix            -       3D connectivity matrix(node*node*subject)
%
% --------------------- OPTIONAL ARGUMENTS ----------------------
%
% regionLabels      -       A cell structure containing names of all nodes in your matrix (at this point its only used for plotting)
% nRand             -       number of randomizations used in random networks for normalization (default = 10)
% prev              -       prevalence for prevalence weighted matrices  (default = 0.75)
% groups            -       1D vector with grouplabels  (default = random)
% nos               -       threshold for the number of streamlines (default = 10)
% percentage        -       threshold for getting the top degree nodes (default = 0.1)
% subjectmask       -       Logical to set which subjects to include (default = all a.k.a. a row of ones)
% PlotLocal         -       Logical to set if we want to plot local metrics (default = 0)
% PlotGlobal        -       Logical to set if we want to plot global metrics (default = 0)
% PlotMatrices      -       Logical to set if we want to plot group average matrices (default = 0)
%
% ---------------------------- OUTPUT ----------------------------
% nResult           -       Graph metrics from binary matrices including all NOS's > nos
% wResult           -       Graph metrics from weighted matrices based on NOS
% pResult           -       Graph metrics from matrices weighted based on group prevalence > prev
% p                 -       Also return all the input settings to check
 
%% check all the inputs and if they do not exist then revert to default settings
% set the larger defaults in case they are not specified
groupdefault = round(rand(1,size(Matrix,3))*3)'; % generate some random numbers
regionsdefault = num2cell([1:size(Matrix,1)])'; % generate a set of numbered labels
subjectsdefault = ones(1,size(Matrix,3));
% input parsing settings
p = inputParser;
p.CaseSensitive = true;
p.Parameters;
p.Results;
p.KeepUnmatched = true;
% set the desired and optional input arguments
addRequired(p,'Matrix',@isnumeric);
addOptional(p,'regionLabels',regionsdefault,@iscell);
addOptional(p,'nRand',10,@isnumeric);
addOptional(p,'PlotLocal',0,@isnumeric);
addOptional(p,'PlotGlobal',0,@isnumeric);
addOptional(p,'PlotMatrices',0,@isnumeric);
addOptional(p,'prev',0.75,@isnumeric);
addOptional(p,'groups',groupdefault, @isnumeric);
addOptional(p,'subjectmask', subjectsdefault, @isnumeric);
addOptional(p,'nos',10, @isnumeric);
addOptional(p,'percentage',0.1, @isnumeric);
% parse the input
parse(p,varargin{:});
% then set/get all the inputs out of this structure
nRand = p.Results.nRand; regionLabels = p.Results.regionLabels; PlotLocal = p.Results.PlotLocal; PlotGlobal = p.Results.PlotGlobal; prev = p.Results.prev; groups = p.Results.groups; Adj = p.Results.Matrix;
nos = p.Results.nos; PlotMatrices = p.Results.PlotMatrices; subjectmask = p.Results.subjectmask; percentage = p.Results.percentage;
 
% first only select those subjects that are included in the subjectmask if
% it is specified (and otherwise the default is to use all)
Adj = Adj(:,:,logical(subjectmask));
groups = groups(logical(subjectmask));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Unweighted/Binary NOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate degree, density, length, and clustering for each subject
% Everything here is done on a single-level, only use prevalence-thresholds
% when looking at group average connectomes
h = waitbar(0,'Running analysis on binary NOS measure'); 
nSubjects = size(Adj,3);
for i = 1:nSubjects
    waitbar(i/nSubjects);
    A = double(Adj(:,:,i) > nos);
    A = A + triu(A,1)';
    % create non-normalized output
    nResult.deg(i,:) = degrees_und(A); % degree
    nResult.dens(i,:) = density_und(A); %density
    nResult.cpl(i,:) = charpath(distance_bin(A)); %characteristic path length
    [M Q] = modularityConsensusFun(A,1,nRand); %optimized clustering coefficient
    nResult.M(i,:) = M; 
    nResult.clustcoeff(i,:) = mean(clustering_coef_bu(A)); 
    nResult.rich(i,:) = rich_club_bu(A,35); %rich-club coefficient
    nResult.trans(i,:) = transitivity_bu(A); %transitivity
    nResult.assor(i,:) = assortativity_bin(A,0); %assortativity
    nResult.effN(i,:) = efficiency_nodal(A); %efficiency
    nResult.part(i,:) = participation_coef(A,M);
    
    % create a random matrix for normalization
    R = randmio_und(A,nRand);
    Crand = mean(clustering_coef_bu(R)); % get random clustering
    Lrand = charpath(distance_bin(R)); % get path length from random matrix
    nResult.Sigma(i,:)=(nResult.clustcoeff(i,:)./Crand)./(nResult.cpl(i,:)./Lrand); % get small world coefficient
    
    for iR = 1:nRand
        R = randmio_und_connected(A, 10);
        RichRand(iR,:) = rich_club_bu(R,35);
        cpl_random(iR,:) = charpath(distance_bin(R));
        clust_random(iR,:) = mean(clustering_coef_bu(R));
        trans_random(iR,:) = transitivity_bu(R);
    end
    RichRand = mean(RichRand,1);
    CplRand = mean(cpl_random,1);
    ClustCoeffRand = mean(clust_random,1);
    TransRand = mean(trans_random,1);
    
    nResult.Norm.rich(i,:) = nResult.rich(i,:)./RichRand;
    nResult.Norm.cpl(i,:) = nResult.cpl(i,:)/CplRand;
    nResult.Norm.clustcoeff(i,:) = nResult.clustcoeff(i,:)/ClustCoeffRand;
    nResult.Norm.trans(i,:) = nResult.trans(i,:)/TransRand;
    
    [nResult.mask(i,:), nResult.net(:,:,i)] = rb_getTop(nResult.deg(i,:),A,percentage);
end
close(h);
 
if PlotGlobal == 1
    figure;
    subplot(2,2,1); boxplot(nResult.dens,groups,'colorgroup',groups);title('density');
    subplot(2,2,2); boxplot(nResult.cpl,groups,'colorgroup',groups);title('characteristic path length');
    subplot(2,2,3); boxplot(nResult.clust,groups,'colorgroup',groups);title('clustering');
    subplot(2,2,4); boxplot(nResult.trans,groups,'colorgroup',groups);title('transitivity');
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Binary Networks based on NOS','HorizontalAlignment','center','VerticalAlignment', 'top');
end
 
if PlotLocal == 1
    figure;
    subplot(4,1,1); bar(mean(nResult.deg(groups == 0,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 0'); title('Degree');
    subplot(4,1,2); bar(mean(nResult.deg(groups == 1,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 1');
    subplot(4,1,3); bar(mean(nResult.deg(groups == 2,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 2');
    subplot(4,1,4); bar(mean(nResult.deg(groups == 3,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 3');
end
 
if PlotMatrices == 1
    figure;
    subplot(2,2,1); imagesc(double(mean(Adj(:,:,(groups == 0)),3))>nos); title('Group 0'); colorbar; colormap(flipud('gray'));
    subplot(2,2,2); imagesc(double(mean(Adj(:,:,(groups == 1)),3))>nos); title('Group 1'); colorbar; colormap(flipud('gray'));
    subplot(2,2,3); imagesc(double(mean(Adj(:,:,(groups == 2)),3))>nos); title('Group 2'); colorbar; colormap(flipud('gray'));
    subplot(2,2,4); imagesc(double(mean(Adj(:,:,(groups == 3)),3))>nos); title('Group 3'); colorbar; colormap(flipud('gray'));
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Weighted NOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Weighted networks
h = waitbar(0,'Running analysis on weighted networks'); 
for i = 1:nSubjects
    waitbar(i/nSubjects);
    A = squeeze(Adj(:,:,i));
    A = A + triu(A,1)';
    % create non-normalized output
    wResult.deg(i,:) = degrees_und(A); % degree
    wResult.dens(i,:) = density_und(A); %density
    wResult.cpl(i,:) = charpath(distance_wei(A)); %characteristic path length
    wResult.trans(i,:) = transitivity_wu(A); %transitivity
    
    [M Q] = modularityConsensusFun(A,1,nRand); %optimized clustering coefficient
    wResult.M(i,:) = M; 
    wResult.clustcoeff(i,:) = mean(clustering_coef_bu(A)); 
    wResult.part(i,:) = participation_coef(A,M);
    
    for iR = 1:nRand
        R = randmio_und_connected(A, 10);% create a random matrix from original
        cpl_random(iR,:) = charpath(distance_wei(R));
        clust_random(iR,:) = mean(clustering_coef_wu(R));
        trans_random(iR,:) = transitivity_bu(R);
    end
    CplRand = mean(cpl_random,1);
    ClustRand = mean(clust_random,1);
    TransRand = mean(trans_random,1);
    
    wResult.Norm.cpl(i,:) = wResult.cpl(i,:)/CplRand;
    wResult.Norm.clustcoeff(i,:) = wResult.clustcoeff(i,:)/ClustRand;
    wResult.Norm.trans(i,:) = wResult.trans(i,:)/TransRand;
    
    [wResult.mask(i,:), wResult.net(:,:,i)] = rb_getTop(wResult.deg(i,:),A,percentage);
end
close(h)
 
if PlotGlobal == 1
    figure;
    subplot(2,2,1); boxplot(wResult.dens,groups,'colorgroup',groups);title('density');
    subplot(2,2,2); boxplot(wResult.cpl,groups,'colorgroup',groups);title('characteristic path length');
    subplot(2,2,3); boxplot(wResult.clust,groups,'colorgroup',groups);title('clustering');
    subplot(2,2,4); boxplot(wResult.trans,groups,'colorgroup',groups);title('transitivity');
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Weighted Networks','HorizontalAlignment','center','VerticalAlignment', 'top');
end
 
if PlotLocal == 1
    figure;
    subplot(4,1,1); bar(mean(wResult.deg(groups == 0,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 0'); title('Degree');
    subplot(4,1,2); bar(mean(wResult.deg(groups == 1,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 1');
    subplot(4,1,3); bar(mean(wResult.deg(groups == 2,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 2');
    subplot(4,1,4); bar(mean(wResult.deg(groups == 3,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 3');
end
 
if PlotMatrices == 1
    figure;
    subplot(2,2,1); imagesc((mean(Adj(:,:,(groups == 0)),3))); title('Group 0'); colorbar;
    subplot(2,2,2); imagesc((mean(Adj(:,:,(groups == 1)),3))); title('Group 1'); colorbar;
    subplot(2,2,3); imagesc((mean(Adj(:,:,(groups == 2)),3))); title('Group 2'); colorbar;
    subplot(2,2,4); imagesc((mean(Adj(:,:,(groups == 3)),3))); title('Group 3'); colorbar;
end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Unweighted by prevalence %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Network with NOS in prevalence percentge of the entire sample
PrevalenceMask = squeeze(double(Adj(:,:,:)>0)); % get all cells that have some NOS
PrevalenceMask = mean(PrevalenceMask,3); % get the mean of those cells
PrevalenceMask = double(PrevalenceMask>prev); % create a mask based on prevalence

h = waitbar(0,'Running analysis on binary prevalence networks'); 
for i = 1:nSubjects
    waitbar(i/nSubjects);
    A = double(Adj(:,:,i) > nos);
    A = A.*PrevalenceMask;
    A = A + triu(A,1)';
    
    % create non-normalized output
    pResult.deg(i,:) = degrees_und(A); % degree
    pResult.dens(i,:) = density_und(A); %density
    pResult.cpl(i,:) = charpath(distance_bin(A)); %characteristic path length
    [M Q] = modularityConsensusFun(A,1,nRand); %optimized clustering coefficient
    pResult.M(i,:) = M; 
    pResult.clustcoeff(i,:) = mean(clustering_coef_bu(A)); 
    pResult.rich(i,:) = rich_club_bu(A,35); %rich-club coefficient
    pResult.trans(i,:) = transitivity_bu(A); %transitivity
    pResult.assor(i,:) = assortativity_bin(A,0); %assortativity
    pResult.effN(i,:) = efficiency_nodal(A); %efficiency
    pResult.part(i,:) = participation_coef(A,M);
    
    % create a random matrix for normalization
    R = randmio_und(A,nRand);
    Crand = mean(clustering_coef_bu(R)); % get random clustering
    Lrand = charpath(distance_bin(R)); % get path length from random matrix
    pResult.Sigma(i,:)=(pResult.clustcoeff(i,:)./Crand)./(pResult.cpl(i,:)./Lrand); % get small world coefficient
    
    for iR = 1:nRand
        R = randmio_und_connected(A, 10);
        RichRand(iR,:) = rich_club_bu(R,35);
        cpl_random(iR,:) = charpath(distance_bin(R));
        clust_random(iR,:) = mean(clustering_coef_bu(R));
        trans_random(iR,:) = transitivity_bu(R);
    end
    RichRand = mean(RichRand,1);
    CplRand = mean(cpl_random,1);
    ClustRand = mean(clust_random,1);
    TransRand = mean(trans_random,1);
    
    pResult.Norm.rich(i,:) = pResult.rich(i,:)./RichRand;
    pResult.Norm.cpl(i,:) = pResult.cpl(i,:)/CplRand;
    pResult.Norm.clustcoeff(i,:) = pResult.clustcoeff(i,:)/ClustRand;
    pResult.Norm.trans(i,:) = pResult.trans(i,:)/TransRand;
    
    [pResult.mask(i,:), pResult.net(:,:,i)] = rb_getTop(pResult.deg(i,:),A,percentage);
end
close(h)
 
if PlotGlobal == 1
    figure;
    subplot(2,2,1); boxplot(pResult.dens,groups,'colorgroup',groups);title('density');
    subplot(2,2,2); boxplot(pResult.cpl,groups,'colorgroup',groups);title('characteristic path length');
    subplot(2,2,3); boxplot(pResult.clust,groups,'colorgroup',groups);title('clustering');
    subplot(2,2,4); boxplot(pResult.trans,groups,'colorgroup',groups);title('transitivity');
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Weighted according to Prevalence','HorizontalAlignment','center','VerticalAlignment', 'top');
end
 
if PlotLocal == 1
    figure;
    subplot(4,1,1); bar(mean(pResult.part(groups == 0,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 0'); title('Degree');
    subplot(4,1,2); bar(mean(pResult.part(groups == 1,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 1');
    subplot(4,1,3); bar(mean(pResult.part(groups == 2,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 2');
    subplot(4,1,4); bar(mean(pResult.part(groups == 3,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 3');
end
 
if PlotMatrices == 1
    figure;
    subplot(2,2,1); imagesc(double(mean(Adj(:,:,(groups == 0)),3).*PrevalenceMask)>0); title('Group 0'); colorbar; colormap(flipud('gray'));
    subplot(2,2,2); imagesc(double(mean(Adj(:,:,(groups == 1)),3).*PrevalenceMask)>0); title('Group 1'); colorbar; colormap(flipud('gray'));
    subplot(2,2,3); imagesc(double(mean(Adj(:,:,(groups == 2)),3).*PrevalenceMask)>0); title('Group 2'); colorbar; colormap(flipud('gray'));
    subplot(2,2,4); imagesc(double(mean(Adj(:,:,(groups == 3)),3).*PrevalenceMask)>0); title('Group 3'); colorbar; colormap(flipud('gray'));
end
 
end
