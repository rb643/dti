function [nResult wResult pResult] = dti(Matrix,varargin)
% 
% For example:
% [nResult wResult pResult] = dti('Matrix',adjacencymatrix);
% [nResult wResult pResult] = dti('Matrix',adjacencymatrix,'regionLabels',regions,'prev',0.75,'PlotLocal',1,'PlotGlobal',1,'nos',10);
% 
% ---------------------------- INPUT ----------------------------
%
% Matrix            -       3D connectivity matrix(node*node*subject)
%
% --------------------- OPTIONAL ARGUMENTS ----------------------
%
% regionLabels      -       A cell structure containing names of all nodes in your matrix (at this point its only used for plotting)
% nRand             -       number of randomizations used in random networks for normalization
% prev              -       prevalence for prevalence weighted matrices
% groups            -       1D vector with grouplabels
% nos               -       threshold for the number of streamlines
% PlotLocal         -       Logical to set if we want to plot local metrics
% PlotGlobal        -       Logical to set if we want to plot global metrics
% PlotMatrices      -       Logical to set if we want to plot group average matrices
%
% ---------------------------- OUTPUT ----------------------------
% nResult           -       Graph metrics from binary matrices including all NOS's > nos
% wResult           -       Graph metrics from weighted matrices based on NOS
% pResult           -       Graph metrics from matrices weighted based on group prevalence > prev

% TODO:
% include subject mask to tell which subject to include (although this is probably better to do beforehand...)

%% check all the inputs and if they do not exist then revert to default settings
% set the larger defaults in case they are not specified
groupdefault = round(rand(1,size(Matrix,3))*3)'; % generate some random numbers
regionsdefault = num2cell([1:size(Matrix,1)])'; % generate a set of numbered labels
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
addOptional(p,'subjectmask', @isnumeric);
addOptional(p,'nos',0, @isnumeric);
% parse the input
parse(p,varargin{:});
% then set get all the inputs out of a structure
nRand = p.Results.nRand; regionLabels = p.Results.regionLabels; PlotLocal = p.Results.PlotLocal; PlotGlobal = p.Results.PlotGlobal; prev = p.Results.prev; groups = p.Results.groups; Adj = p.Results.Matrix;
nos = p.Results.nos; PlotMatrices = p.Results.PlotMatrices;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Unweighted NOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate degree, density, length, and clustering for each subject
% Everything here is done on a single-level, only use prevalence-thresholds
% when looking at group average connectomes
nSubjects = size(Adj,3);
for i = 1:nSubjects
    A = double(Adj(:,:,i) > nos);
    A = A + triu(A,1)';
    % create non-normalized output
    nResult.deg(i,:) = degrees_und(A); % degree
    nResult.dens(i,:) = density_und(A); %density
    nResult.cpl(i,:) = charpath(distance_bin(A)); %characteristic path length
    nResult.clust(i,:) = mean(clustering_coef_bu(A)); %clustering coefficient
    nResult.rich(i,:) = rich_club_bu(A,35); %rich-club coefficient
    nResult.trans(i,:) = transitivity_bu(A); %transitivity
    nResult.assor(i,:) = assortativity_bin(A,0); %assortativity
    nResult.effN(i,:) = efficiency_nodal(A); %efficiency
    
    % create a random matrix for normalization
    R = randmio_und(A,nRand);
    Crand = mean(clustering_coef_bu(R)); % get random clustering
    Lrand = charpath(distance_bin(R)); % get path length from random matrix
    nResult.Sigma(i,:)=(nResult.clust(i,:)./Crand)./(nResult.cpl(i,:)./Lrand); % get small world coefficient
    cpl_random(i,:) = mean(squareform(distance_bin(R)));
    clust_random(i,:) = mean(clustering_coef_bu(R));
    for iR = 1:nRand
        TempClub(iR,:,:) = randmio_und_connected(A, nRand);% create a random matrix from original
        RichRand(iR,:) = rich_club_bu(squeeze(TempClub(iR,:,:)),35);
    end
    RichRand = mean(RichRand,1);
    
    nResult.Norm.rich(i,:) = nResult.rich(i,:)./RichRand;
end

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
for i = 1:nSubjects
    A = squeeze(Adj(:,:,i));
    A = A + triu(A,1)';
    % create non-normalized output
    wResult.deg(i,:) = degrees_und(A); % degree
    wResult.dens(i,:) = density_und(A); %density
    wResult.cpl(i,:) = charpath(distance_wei(A)); %characteristic path length
    wResult.clust(i,:) = mean(clustering_coef_wu(A)); %clustering coefficient
    wResult.trans(i,:) = transitivity_wu(A); %transitivity
end

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

for i = 1:nSubjects
    A = double(Adj(:,:,i) > nos);
    A = A.*PrevalenceMask;
    A = A + triu(A,1)';
  
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
    subplot(4,1,1); bar(mean(pResult.deg(groups == 0,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 0'); title('Degree');
    subplot(4,1,2); bar(mean(pResult.deg(groups == 1,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 1');
    subplot(4,1,3); bar(mean(pResult.deg(groups == 2,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 2');
    subplot(4,1,4); bar(mean(pResult.deg(groups == 3,:)),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ylabel('Group 3');
end

if PlotMatrices == 1
    figure;
    subplot(2,2,1); imagesc(double(mean(Adj(:,:,(groups == 0)),3).*PrevalenceMask)>0); title('Group 0'); colorbar; colormap(flipud('gray'));
    subplot(2,2,2); imagesc(double(mean(Adj(:,:,(groups == 1)),3).*PrevalenceMask)>0); title('Group 1'); colorbar; colormap(flipud('gray'));
    subplot(2,2,3); imagesc(double(mean(Adj(:,:,(groups == 2)),3).*PrevalenceMask)>0); title('Group 2'); colorbar; colormap(flipud('gray'));
    subplot(2,2,4); imagesc(double(mean(Adj(:,:,(groups == 3)),3).*PrevalenceMask)>0); title('Group 3'); colorbar; colormap(flipud('gray'));
end

end
