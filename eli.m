function [Result wResult pResult] = eli(Matrix, nRand, PlotAverage, PlotAll, GroupMask, prev, metric)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matrix            -       4D connectivity matrix(node*node*measure*subject)
%% UseMask           -       1D logical with which subjects to include, if not supplied defaults to all
%% nRand             -       number of randomization used in random networks for normalization
%% PlotAverage       -       Logical to set if we want to plot average
%% PlotAll           -       Logical to set if we want to plot the output of some graph metrics
%% GroupMask         -       1D vector with grouplabels
%% prev              -       prevalence for prevalence weighted matrices
%% metric            -       metric that we want to use to construct matrices (defaults to number of streamlines)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Result            -       Graph metrics from binary matrices including all NOS's >0
%% wResult           -       Graph metrics from weighted matrices 
%% pResult           -       Graph metrics from matrices weighted based on group prevalence

if exist('GroupMask')
    % Assign default group membership
    groups = GroupMask;
else
    groups = [1 2 3 0	2 0	0	0	3	0	2	3	0	1	1	0	0	0	1	0	2	2	2	1	2	3	3	3	3	2	1	3	2	3	0	0	2	0	1	2	0	0	0	1	1	1	3	2	3	2	1	3	1	1	1	2	1	1	2	1	1	2	1	1	1	2	1	2	1	0	2	1	3	1	1	0	0	3	0	0	0	0	2	1	2	3	2	2	1	1	6	1	3	3	2	3	3]';
end

if exist('nRand'); nRand = nRand; else nRand =10; end
if exist('PlotAverage'); PlotAverage = PlotAverage; else PlotAverage = 0; end
if exist('PlotAverage'); PlotAll = PlotAll; else PlotAll = 0; end
if exist('metric'); metric = metric; else metric = 1; end

% Remove subjects with missing clinical data: CAT118, CAT120, CAT122
connectivity = Matrix;
nSubjects = size(connectivity,4);
% Extract NOS matrix for every subject and remove singleton dimensions
NOS_w = squeeze(connectivity(:,:,metric,:));

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
P_Sample = mean(connectivity(:,:,metric,:) > 0, 4);
A_Sample = double(P>=0.6);

P_0= mean(connectivity(:,:,metric,Group0) > 0, 4);
A_0 = double(P_0>=0.5);

P_1 = mean(connectivity(:,:,metric,Group1) > 0, 4);
A_1 = double(P_1>=0.5);

P_2= mean(connectivity(:,:,metric,Group2) > 0, 4);
A_2 = double(P_2>=0.5);

P_3 = mean(connectivity(:,:,metric,Group3) > 0, 4);
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
end

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

end

