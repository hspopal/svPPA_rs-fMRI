%% Add relevent paths

%addpath('/Users/tuk12127/matlab/barwitherr')
%addpath('/Users/tuk12127/matlab/xticklabel_rotate')
%addpath('/Users/haroonpopal/matlab/computeCohen_d')

%% Import ROI to ROI correltions for each subject
cd '/Users/haroonpopal/Google_Drive/svPPA_rs-fMRI/results/matlab'

% Create structures for each group
% Clear after each script is run so that the groups do not get mixed
% run('indiv_subj_corrs_z_YC.m')
% clear
% run('indiv_subj_corrs_z_OC.m')
% clear
% run('indiv_subj_corrs_z_svPPA_n12.m')
% clear

load('indiv_subj_corrs_z_YC.mat')
load('indiv_subj_corrs_z_OC.mat')
load('indiv_subj_corrs_z_svPPA_n12.mat')



%% Create an array with group identifiers in each cell
YC_length = length(lMT_corrs_YC.lFEF_lMT);
OC_length = length(lMT_corrs_OC.lFEF_lMT);
svPPA_length = length(lMT_corrs_svPPA.lFEF_lMT);

Groups_array = cell(1,YC_length+OC_length+svPPA_length);
Groups_array(1:YC_length) = {'YC'};
Groups_array(YC_length+1:YC_length+OC_length) = {'OC'};
Groups_array(YC_length+OC_length+1:YC_length+OC_length+svPPA_length) = {'svPPA'};

% These ROIs need to be named exactly as they were in the config file used
% to run fc_analysis in order to create the correct corrs_all structure
% belows. Afterwards, the ROI names can be changed for figures, etc.
ROIs = {'lFEF','lSPL','lMT','lV1','lV1c','lV1p','lOFA','lFFA','lATC','lmPFC','lMTG','lAG','lPCC','rFEF','rSPL','rMT','rV1','rV1c','rV1p','rOFA','rFFA','rATC','rMTG','rmPFC','rAG','rPCC'};

%% create _corrs_all for all ROIs

all_corrs = struct;
for j=1:length(ROIs)
    strc_YC = eval(char(strcat(ROIs(j), '_corrs_YC')));
    strc_OC = eval(char(strcat(ROIs(j), '_corrs_OC')));
    strc_svPPA = eval(char(strcat(ROIs(j), '_corrs_svPPA')));
    SNames = fieldnames(strc_YC);
    for i=1:(length(ROIs)-1);
        all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(SNames{i}) = [strc_YC.(SNames{i}); strc_OC.(SNames{i}); strc_svPPA.(SNames{i})];
    end
end


%% One-way ANOVAs

% structure with means, SE, and p-values for each ROI, by column (ie.
% column 1 is means, 2 is SE, 3 is p-values
% p-values are output from anova1 (with rows being 1v2, 1v3, 2v3)
for j=1:length(ROIs)
    for i=1:length(ROIs)
        if ismember(ROIs(j),ROIs(i))
            continue
        else
            if isfield(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))), strcat(ROIs(j), '_', ROIs(i)))
                temp_pair = strcat(ROIs(j), '_', ROIs(i));
                [p,tbl,stats] = anova1(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(char(temp_pair)), Groups_array);
                [c,m] = multcompare(stats);
                corr_anova1_analyses.(char(ROIs(j))).(char(ROIs(i))) = [m, [c(1,6); c(2,6); c(3,6)]];
                close all
            elseif isfield(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))), strcat(ROIs(i), '_', ROIs(j)))
                temp_pair = strcat(ROIs(i), '_', ROIs(j));
                [p,tbl,stats] = anova1(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(char(temp_pair)), Groups_array);
                [c,m] = multcompare(stats);
                corr_anova1_analyses.(char(ROIs(j))).(char(ROIs(i))) = [m, [c(1,6); c(2,6); c(3,6)]];
                close all
            else
                disp(strcat(ROIs(i), ' is not present in  ', ROIs(j)))
            end
        end
    end
end


%% Create structures for means, SEs, and ROI names, to easily create bar graphs

corr_means = struct;
corr_SEs = struct;
corr_ROIs = struct;
for j=1:length(ROIs)
    for i=1:length(ROIs)
        if ismember(ROIs(j),ROIs(i))
            continue
        else
            if isfield(corr_means, (char(ROIs(j))))
                corr_means.(char(ROIs(j))) = [corr_means.(char(ROIs(j))); corr_anova1_analyses.(char(ROIs(j))).(char(ROIs(i)))(:,1)'];
                corr_SEs.(char(ROIs(j))) = [corr_SEs.(char(ROIs(j))); corr_anova1_analyses.(char(ROIs(j))).(char(ROIs(i)))(:,2)'];
                corr_ROIs.(char(ROIs(j))) = [corr_ROIs.(char(ROIs(j))); ROIs(i)];
            else
                corr_means.(char(ROIs(j))) = [];
                corr_means.(char(ROIs(j))) = [corr_anova1_analyses.(char(ROIs(j))).(char(ROIs(i)))(:,1)'];
                corr_SEs.(char(ROIs(j))) = [];
                corr_SEs.(char(ROIs(j))) = [corr_anova1_analyses.(char(ROIs(j))).(char(ROIs(i)))(:,2)'];
                corr_ROIs.(char(ROIs(j))) = [];
                corr_ROIs.(char(ROIs(j))) = ROIs(i);
            end
        end
    end
end


%% Create bar graphs

% inter-network RH
cd '/autofs/cluster/animal/users/hsp13/projects/svPPA_rs-fMRI/results/matlab/DAN_VS_DMN/figures/bar_graphs/inter_network_5mm'
for i=12:length(ROIs);
    h = figure(i);
    b_err = barwitherr(corr_SEs.(char(ROIs(i)))(12:end,:), corr_means.(char(ROIs(i)))(12:end,:));
    b_err(1).FaceColor = 'b';
    b_err(2).FaceColor = 'r';
    b_err(3).FaceColor = 'g';
    legend('YC','OC','svPPA','Location', 'eastoutside')
    xticklabel_rotate([1:11],45,corr_ROIs.(char(ROIs(i)))(12:end))
    ylabel('Mean z-score');
    title(strcat('Mean  ', ROIs(i), ' Correlations'));
    saveas(h,strcat(ROIs{i},'_zscores.png'));
end
close all

% inter-network LH
cd '/autofs/cluster/animal/users/hsp13/projects/svPPA_rs-fMRI/results/matlab/DAN_VS_DMN/figures/bar_graphs/inter_network_5mm'
for i=1:11;
    h = figure(i);
    b_err = barwitherr(corr_SEs.(char(ROIs(i)))([1:10 22],:), corr_means.(char(ROIs(i)))([1:10 22],:));
    b_err(1).FaceColor = 'b';
    b_err(2).FaceColor = 'r';
    b_err(3).FaceColor = 'g';
    legend('YC','OC','svPPA','Location', 'eastoutside')
    xticklabel_rotate([1:11],45,corr_ROIs.(char(ROIs(i)))([1:10 22]))
    ylabel('Mean z-score');
    title(strcat('Mean  ', ROIs(i), ' Correlations'));
    saveas(h,strcat(ROIs{i},'_zscores.png'));
end
close all

cd '/autofs/cluster/animal/users/hsp13/projects/svPPA_rs-fMRI/results/matlab/DAN_VS_DMN/figures/bar_graphs/inter_network_5mm'
figure;
b_err = barwitherr([corr_SEs.PCC(1:11,:)], [corr_means.PCC(1:11,:)]);
b_err(1).FaceColor = 'b';
b_err(2).FaceColor = 'r';
b_err(3).FaceColor = 'g';
legend('YC','OC','SD','Location', 'eastoutside')
xticklabel_rotate([1:11],45,[corr_ROIs.PCC(1:11)])
ylabel('Mean z-score', 'FontSize', 12);
%ylim([-.1 .45])
title('Mean PCC Correlations');
saveas(h,'lPCC_zscores.png');

%% Create bar graphs OC vs svPPA

% inter-network RH
cd '/Users/tuk12127/Google_Drive/svPPA_rs-fMRI/results/matlab/DAN_VS/figures/bar_graphs/inter_network_5mm'
for i=8:length(ROIs);
    h = figure(i);
    b_err = barwitherr(corr_SEs.(char(ROIs(i)))(8:end,2:3), corr_means.(char(ROIs(i)))(8:end,2:3));
    b_err(1).FaceColor = 'b';
    b_err(2).FaceColor = 'r';
    legend('OC','svPPA','Location', 'eastoutside')
    xticklabel_rotate([1:6],45,corr_ROIs.(char(ROIs(i)))(8:end))
    ylabel('Mean Correlation (r)');
    title(strcat('Mean  ', ROIs(i), ' Correlations'));
    saveas(h,strcat(ROIs{i},'_OC_svPPA_zscores.png'));
end
close all


b_err = barwitherr(corr_SEs.rMTG(12:end,2:3), corr_means.rMTG(12:end,2:3));
b_err(1).FaceColor = 'b';
b_err(2).FaceColor = 'r';
legend('OC','svPPA','Location', 'eastoutside')
xticklabels(corr_ROIs.rMTG(12:end))
ylabel('Mean Correlation (z)');
title(strcat('Mean rMTG Correlations'));


%% Conduct one-tailed and two-tailed t-tests
%   variances of groups are mostly equal
%   one-tailed should be used for DAN to VisAN ROIs
%   two-tailed should be used for DAN to DMN ROIs (not clear hypotheses)
for j=12:length(ROIs)
    for i=12:length(ROIs)
        if ismember(ROIs(j),ROIs(i))
            continue
        else
            if isfield(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))), strcat(ROIs(j), '_', ROIs(i)))
                temp_pair = strcat(ROIs(j), '_', ROIs(i));
                [h, p, ci, stats] = ttest2(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(char(temp_pair))(find(strcmp(Groups_array, 'OC'))), all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(char(temp_pair))(find(strcmp(Groups_array, 'svPPA'))),'Tail','Right');
                corr_ttest_analyses.(char(ROIs(j))).(char(ROIs(i))) = [h, p, ci(2), stats.tstat];
                [h, p, ci, stats] = ttest2(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(char(temp_pair))(find(strcmp(Groups_array, 'OC'))), all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(char(temp_pair))(find(strcmp(Groups_array, 'svPPA'))));
                corr_ttest_2tail_analyses.(char(ROIs(j))).(char(ROIs(i))) = [h, p, ci(2), stats.tstat];
                close all
            elseif isfield(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))), strcat(ROIs(i), '_', ROIs(j)))
                temp_pair = strcat(ROIs(i), '_', ROIs(j));
                [h, p, ci, stats] = ttest2(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(char(temp_pair))(find(strcmp(Groups_array, 'OC'))), all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(char(temp_pair))(find(strcmp(Groups_array, 'svPPA'))),'Tail','Right');
                corr_ttest_analyses.(char(ROIs(j))).(char(ROIs(i))) = [h, p, ci(2), stats.tstat];
                [h, p, ci, stats] = ttest2(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(char(temp_pair))(find(strcmp(Groups_array, 'OC'))), all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(char(temp_pair))(find(strcmp(Groups_array, 'svPPA'))));
                corr_ttest_2tail_analyses.(char(ROIs(j))).(char(ROIs(i))) = [h, p, ci(2), stats.tstat];
                close all
            else
                disp(strcat(ROIs(i), ' is not present in rSPL'))
            end
        end
    end
end

%% ROI t-tests and effect sizes
[h, p, ci, stats] = ttest2(rSPL_corrs_OC.rFFA_rSPL, rSPL_corrs_svPPA.rFFA_rSPL);

% Cohen's d
n1 = length(rSPL_corrs_OC.rFFA_rSPL);
n2 = length(rSPL_corrs_svPPA.rFFA_rSPL);
var1 = std(rSPL_corrs_OC.rFFA_rSPL)^2;
var2 = std(rSPL_corrs_svPPA.rFFA_rSPL)^2;
Spooled = sqrt(((n1-1)*var1 + (n2-1)*var2) / ((n1-1) + (n2-1)));
cohen_d = (mean(rSPL_corrs_OC.rFFA_rSPL) - mean(rSPL_corrs_svPPA.rFFA_rSPL)) / Spooled;

cohend.rSPL.rFFA = computeCohen_d(rSPL_corrs_OC.rFFA_rSPL,rSPL_corrs_svPPA.rFFA_rSPL,'independent');
cohend.rSPL.rV1  = computeCohen_d(rSPL_corrs_OC.rV1_rSPL,rSPL_corrs_svPPA.rV1_rSPL,'independent');
cohend.rSPL.rAG  = computeCohen_d(rSPL_corrs_OC.rAG_rSPL,rSPL_corrs_svPPA.rAG_rSPL,'independent');
cohend.rSPL.rPCC = computeCohen_d(rSPL_corrs_OC.rPCC_rSPL,rSPL_corrs_svPPA.rPCC_rSPL,'independent');
cohend.rMTG.rSPL = computeCohen_d(rMTG_corrs_OC.rMTG_rSPL,rMTG_corrs_svPPA.rMTG_rSPL,'independent');
cohend.rMTG.rMT  = computeCohen_d(rMTG_corrs_OC.rMTG_rMT,rMTG_corrs_svPPA.rMTG_rMT,'independent');
cohend.rMTG.rAG  = computeCohen_d(rMTG_corrs_OC.rMTG_rAG,rMTG_corrs_svPPA.rMTG_rAG,'independent');
cohend.rMTG.rPCC  = computeCohen_d(rMTG_corrs_OC.rPCC_rMTG,rMTG_corrs_svPPA.rPCC_rMTG,'independent');
cohend.rFEF.rFFA  = computeCohen_d(rFEF_corrs_OC.rFFA_rFEF,rFEF_corrs_svPPA.rFFA_rFEF,'independent');
cohend.rMT.rATC  = computeCohen_d(rMT_corrs_OC.rATC_rMT,rMT_corrs_svPPA.rATC_rMT,'independent');


%% One-sample t-tests


ROIs = {'lFEF','lSPL','lMT','lV1c','lOFA','lFFA','lATC','lmPFC','lMTG','lAG','lPCC','rFEF','rSPL','rMT','rV1c','rOFA','rFFA','rATC','rMTG','rmPFC','rAG','rPCC'};


for j=1:length(ROIs)
    for i=1:length(ROIs)
        if ismember(ROIs(j),ROIs(i))
            continue
        else
            if isfield(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))), strcat(ROIs(j), '_', ROIs(i)))
                temp_pair = strcat(ROIs(j), '_', ROIs(i));
                [h, p, ci, stats] = ttest(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(char(temp_pair))(find(strcmp(Groups_array, 'YC'))));
                corr_ttest_1samp_analyses_YC.(char(ROIs(j))).(char(ROIs(i))) = [h, p, ci(2), stats.tstat];
                [h, p, ci, stats] = ttest(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(char(temp_pair))(find(strcmp(Groups_array, 'OC'))));
                corr_ttest_1samp_analyses_OC.(char(ROIs(j))).(char(ROIs(i))) = [h, p, ci(2), stats.tstat];
                [h, p, ci, stats] = ttest(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(char(temp_pair))(find(strcmp(Groups_array, 'svPPA'))));
                corr_ttest_1samp_analyses_svPPA.(char(ROIs(j))).(char(ROIs(i))) = [h, p, ci(2), stats.tstat];
                close all
            elseif isfield(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))), strcat(ROIs(i), '_', ROIs(j)))
                temp_pair = strcat(ROIs(i), '_', ROIs(j));
                [h, p, ci, stats] = ttest(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(char(temp_pair))(find(strcmp(Groups_array, 'YC'))));
                corr_ttest_1samp_analyses_YC.(char(ROIs(j))).(char(ROIs(i))) = [h, p, ci(2), stats.tstat];
                [h, p, ci, stats] = ttest(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(char(temp_pair))(find(strcmp(Groups_array, 'OC'))));
                corr_ttest_1samp_analyses_OC.(char(ROIs(j))).(char(ROIs(i))) = [h, p, ci(2), stats.tstat];
                [h, p, ci, stats] = ttest(all_corrs.(char(strcat(ROIs(j), '_corrs_all'))).(char(temp_pair))(find(strcmp(Groups_array, 'svPPA'))));
                corr_ttest_1samp_analyses_svPPA.(char(ROIs(j))).(char(ROIs(i))) = [h, p, ci(2), stats.tstat];
                close all
            else
                disp(strcat(ROIs(i), ' is not present in rSPL'))
            end
        end
    end
end


temp_nodes = rFFA_corrs_YC.rATC_rFFA;
[h,p,ci,stats] = ttest(temp_nodes)
mean(temp_nodes)

%% Create adjacency Matrix

ROIs = {'lFEF','lSPL','lMT','lV1','lV1c','lV1p','lOFA','lFFA','lATC','lmPFC','lMTG','lAG','lPCC','rFEF','rSPL','rMT','rV1','rV1c','rV1p','rOFA','rFFA','rATC','rMTG','rmPFC','rAG','rPCC'};

% I made the diagonals 0 to threshold so that matrices are not skewed
adj_YC = [0; corr_means.(char(ROIs(1)))(:,1)];
for i=2:length(ROIs)
    adj_YC = [adj_YC [corr_means.(char(ROIs(i)))(1:(i-1),1); 0; corr_means.(char(ROIs(i)))(i:end,1)]];
end

adj_OC = [0; corr_means.(char(ROIs(1)))(:,2)];
for i=2:length(ROIs)
    adj_OC = [adj_OC [corr_means.(char(ROIs(i)))(1:(i-1),2); 0; corr_means.(char(ROIs(i)))(i:end,2)]];
end

adj_svPPA = [0; corr_means.(char(ROIs(1)))(:,3)];
for i=2:length(ROIs)
    adj_svPPA = [adj_svPPA [corr_means.(char(ROIs(i)))(1:(i-1),3); 0; corr_means.(char(ROIs(i)))(i:end,3)]];
end


% Matrices for within hemisphere ROIs
adj_YC_rh = adj_YC(14:end,14:end);
adj_YC_lh = adj_YC(1:13,1:13);

adj_OC_rh = adj_OC(14:end,14:end);
adj_OC_lh = adj_OC(1:13,1:13);

adj_svPPA_rh = adj_svPPA(14:end,14:end);
adj_svPPA_lh = adj_svPPA(1:13,1:13);


% Threshold based on either positive or anti correlations
thresh = 0.1;
adj_YC_thresh = adj_YC_rh;
adj_OC_thresh = adj_OC_rh;
adj_svPPA_thresh = adj_svPPA_rh;
if thresh > 0  % if thresh is positive, threshold everything below
    adj_YC_thresh((adj_YC_thresh < thresh)) = 0;
    adj_OC_thresh((adj_OC_thresh < thresh)) = 0;
    adj_svPPA_thresh((adj_svPPA_thresh < thresh)) = 0;
elseif thresh < 0  % if thresh is negative, threshold everything above
    adj_YC_thresh((adj_YC_thresh > thresh)) = 0;
    adj_OC_thresh((adj_OC_thresh > thresh)) = 0;
    adj_svPPA_thresh((adj_svPPA_thresh > thresh)) = 0;
end

% Remove unwanted rows and columns
adj_YC_thresh(:,6) = [];
adj_YC_thresh(:,4) = [];
adj_YC_thresh(6,:) = [];
adj_YC_thresh(4,:) = [];

adj_OC_thresh(:,6) = [];
adj_OC_thresh(:,4) = [];
adj_OC_thresh(6,:) = [];
adj_OC_thresh(4,:) = [];

adj_svPPA_thresh(:,6) = [];
adj_svPPA_thresh(:,4) = [];
adj_svPPA_thresh(6,:) = [];
adj_svPPA_thresh(4,:) = [];

ROIs([4 6 17 19]) = [];

%% Create adjacency Matrix

%ROIs = {'lFEF','lSPL','lMT','lV1','lV1c','lV1p','lOFA','lFFA','lATC','lmPFC','lMTG','lAG','lPCC','rFEF','rSPL','rMT','rV1','rV1c','rV1p','rOFA','rFFA','rATC','rMTG','rmPFC','rAG','rPCC'};

% I made the diagonals 0 to threshold so that matrices are not skewed
adj_YC = [0; corr_means.(char(ROIs(1)))(:,1)];
for i=2:length(ROIs)
    adj_YC = [adj_YC [corr_means.(char(ROIs(i)))(1:(i-1),1); 0; corr_means.(char(ROIs(i)))(i:end,1)]];
end

adj_OC = [0; corr_means.(char(ROIs(1)))(:,2)];
for i=2:length(ROIs)
    adj_OC = [adj_OC [corr_means.(char(ROIs(i)))(1:(i-1),2); 0; corr_means.(char(ROIs(i)))(i:end,2)]];
end

adj_svPPA = [0; corr_means.(char(ROIs(1)))(:,3)];
for i=2:length(ROIs)
    adj_svPPA = [adj_svPPA [corr_means.(char(ROIs(i)))(1:(i-1),3); 0; corr_means.(char(ROIs(i)))(i:end,3)]];
end


% Matrices for within hemisphere ROIs
adj_YC_rh = adj_YC(14:end,14:end);
adj_YC_lh = adj_YC(1:13,1:13);

adj_OC_rh = adj_OC(14:end,14:end);
adj_OC_lh = adj_OC(1:13,1:13);

adj_svPPA_rh = adj_svPPA(14:end,14:end);
adj_svPPA_lh = adj_svPPA(1:13,1:13);







3%% Graph Theory

% Use this ROI list to change the name of ROIs in the figures
%ROIs = {'lFEF','lSPL','lMT','lV1c','lOFA','lFFA','lPRC','lmPFC','lMTG','lAG','lPCC','rFEF','rSPL','rMT','rV1c','rOFA','rFFA','rPRC','rMTG','rmPFC','rAG','rPCC'};
ROIs(7) = {'lPRC'};
ROIs(18) = {'rPRC'};


% Biograph - to create simple graphs in matlab
graph_YC = biograph(adj_YC_thresh,ROIs(12:end),'EdgeType','straight','ShowArrows','off','LayoutType','equilibrium','NodeAutoSize','off');
graph_OC = biograph(adj_OC_thresh,ROIs(12:end),'EdgeType','straight','ShowArrows','off','LayoutType','equilibrium','NodeAutoSize','off');
graph_svPPA = biograph(adj_svPPA_thresh,ROIs(12:end),'EdgeType','straight','ShowArrows','off','LayoutType','equilibrium','NodeAutoSize','off');

% Change the graph line sizes for easier visualization and change the ROI
% colors to represent networks
% Colors are defined by the RGB scale from 0 to 1
for i=1:length(graph_YC.Edges),
    graph_YC.Edges(i).LineWidth = abs(graph_YC.Edges(i).Weight*10);
    if i < 4,
        graph_YC.Nodes(i).Color = [0,1,0];
        graph_YC.Nodes(i).LineColor = [0,1,0];
    elseif i == 4,
        graph_YC.Nodes(i).Color = [1,1,1];
        graph_YC.Nodes(i).LineColor = [0,0,0];
    elseif i < 8,
        graph_YC.Nodes(i).Color = [1,0,1];
        graph_YC.Nodes(i).LineColor = [1,0,1];
    elseif i < 12
        graph_YC.Nodes(i).Color = [1,0.4196,0.0078];
        graph_YC.Nodes(i).LineColor = [1,0.4196,0.0078];
    end
    if i < 12,
        graph_YC.Nodes(i).Size = [30 30];
        graph_YC.Nodes(i).Shape = 'circle';
        graph_YC.Nodes(i).FontSize = 20;
    end
end
for i=1:length(graph_OC.Edges),
    graph_OC.Edges(i).LineWidth = abs(graph_OC.Edges(i).Weight*10);
    if i < 4,
        graph_OC.Nodes(i).Color = [0,1,0];
        graph_OC.Nodes(i).LineColor = [0,1,0];
    elseif i == 4,
        graph_OC.Nodes(i).Color = [1,1,1];
        graph_OC.Nodes(i).LineColor = [0,0,0];
    elseif i < 8,
        graph_OC.Nodes(i).Color = [1,0,1];
        graph_OC.Nodes(i).LineColor = [1,0,1];
    elseif i < 12
        graph_OC.Nodes(i).Color = [1,0.4196,0.0078];
        graph_OC.Nodes(i).LineColor = [1,0.4196,0.0078];
    end
    if i < 12,
        graph_OC.Nodes(i).Size = [30 30];
        graph_OC.Nodes(i).Shape = 'circle';
        graph_OC.Nodes(i).FontSize = 20;
    end
end
for i=1:length(graph_svPPA.Edges),
    graph_svPPA.Edges(i).LineWidth = abs(graph_svPPA.Edges(i).Weight*10);
    if i < 4,
        graph_svPPA.Nodes(i).Color = [0,1,0];
        graph_svPPA.Nodes(i).LineColor = [0,1,0];
    elseif i == 4,
        graph_svPPA.Nodes(i).Color = [1,1,1];
        graph_svPPA.Nodes(i).LineColor = [0,0,0];
    elseif i < 8,
        graph_svPPA.Nodes(i).Color = [1,0,1];
        graph_svPPA.Nodes(i).LineColor = [1,0,1];
    elseif i < 12
        graph_svPPA.Nodes(i).Color = [1,0.4196,0.0078];
        graph_svPPA.Nodes(i).LineColor = [1,0.4196,0.0078];
    end
    if i < 12,
        graph_svPPA.Nodes(i).Size = [30 30];
        graph_svPPA.Nodes(i).Shape = 'circle';
        graph_svPPA.Nodes(i).FontSize = 20;
    end
end

view(graph_YC)
view(graph_OC)
view(graph_svPPA)


%% Graph thesholded by significant t-tests

ROIs = {'rFEF','rSPL','rMT','rV1c','rOFA','rFFA','rATC','rMTG','rmPFC','rAG','rPCC'};


% I made the diagonals 0 to threshold so that matrices are not skewed
adj_YC_ttest = adj_YC_thresh;
adj_OC_ttest = adj_OC_thresh;
adj_svPPA_ttest = adj_svPPA_thresh;
for j=1:length(ROIs)
    for i=1:length(ROIs)
        if ismember(ROIs(j),ROIs(i))
            continue
        end
        if corr_ttest_1samp_analyses_YC.(char(ROIs(j))).(char(ROIs(i)))(2) < 0.00091
            continue
        else
            adj_YC_ttest(j,i) = 0;
        end
        if corr_ttest_1samp_analyses_OC.(char(ROIs(j))).(char(ROIs(i)))(2) < 0.00091
            continue
        else
            adj_OC_ttest(j,i) = 0;
        end
        if corr_ttest_1samp_analyses_svPPA.(char(ROIs(j))).(char(ROIs(i)))(2) < 0.00091
            continue
        else
            adj_svPPA_ttest(j,i) = 0;
        end
    end
end


% Biograph - to create simple graphs in matlab
graph_YC = biograph(adj_YC_ttest,ROIs,'EdgeType','straight','ShowArrows','off','LayoutType','equilibrium','NodeAutoSize','off');
graph_OC = biograph(adj_OC_ttest,ROIs,'EdgeType','straight','ShowArrows','off','LayoutType','equilibrium','NodeAutoSize','off');
graph_svPPA = biograph(adj_svPPA_ttest,ROIs,'EdgeType','straight','ShowArrows','off','LayoutType','equilibrium','NodeAutoSize','off');


% Change the graph line sizes for easier visualization and change the ROI
% colors to represent networks
% Colors are defined by the RGB scale from 0 to 1
for i=1:length(graph_YC.Edges),
    graph_YC.Edges(i).LineWidth = abs(graph_YC.Edges(i).Weight*10);
    if i < 4,
        graph_YC.Nodes(i).Color = [0,1,0];
        graph_YC.Nodes(i).LineColor = [0,1,0];
    elseif i == 4,
        graph_YC.Nodes(i).Color = [1,1,1];
        graph_YC.Nodes(i).LineColor = [0,0,0];
    elseif i < 8,
        graph_YC.Nodes(i).Color = [1,0,1];
        graph_YC.Nodes(i).LineColor = [1,0,1];
    elseif i < 12
        graph_YC.Nodes(i).Color = [1,0.4196,0.0078];
        graph_YC.Nodes(i).LineColor = [1,0.4196,0.0078];
    end
    if i < 12,
        graph_YC.Nodes(i).Size = [30 30];
        graph_YC.Nodes(i).Shape = 'circle';
        graph_YC.Nodes(i).FontSize = 20;
    end
end
for i=1:length(graph_OC.Edges),
    graph_OC.Edges(i).LineWidth = abs(graph_OC.Edges(i).Weight*10);
    if i < 4,
        graph_OC.Nodes(i).Color = [0,1,0];
        graph_OC.Nodes(i).LineColor = [0,1,0];
    elseif i == 4,
        graph_OC.Nodes(i).Color = [1,1,1];
        graph_OC.Nodes(i).LineColor = [0,0,0];
    elseif i < 8,
        graph_OC.Nodes(i).Color = [1,0,1];
        graph_OC.Nodes(i).LineColor = [1,0,1];
    elseif i < 12
        graph_OC.Nodes(i).Color = [1,0.4196,0.0078];
        graph_OC.Nodes(i).LineColor = [1,0.4196,0.0078];
    end
    if i < 12,
        graph_OC.Nodes(i).Size = [30 30];
        graph_OC.Nodes(i).Shape = 'circle';
        graph_OC.Nodes(i).FontSize = 20;
    end
end
for i=1:length(graph_svPPA.Edges),
    graph_svPPA.Edges(i).LineWidth = abs(graph_svPPA.Edges(i).Weight*10);
    if i < 4,
        graph_svPPA.Nodes(i).Color = [0,1,0];
        graph_svPPA.Nodes(i).LineColor = [0,1,0];
    elseif i == 4,
        graph_svPPA.Nodes(i).Color = [1,1,1];
        graph_svPPA.Nodes(i).LineColor = [0,0,0];
    elseif i < 8,
        graph_svPPA.Nodes(i).Color = [1,0,1];
        graph_svPPA.Nodes(i).LineColor = [1,0,1];
    elseif i < 12
        graph_svPPA.Nodes(i).Color = [1,0.4196,0.0078];
        graph_svPPA.Nodes(i).LineColor = [1,0.4196,0.0078];
    end
    if i < 12,
        graph_svPPA.Nodes(i).Size = [30 30];
        graph_svPPA.Nodes(i).Shape = 'circle';
        graph_svPPA.Nodes(i).FontSize = 20;
    end
end

view(graph_YC)
view(graph_OC)
view(graph_svPPA)



%% Adjacency Matrix Plots

% Create adjacency Matrix

ROIs = {'lFEF','lSPL','lMT','lV1','lV1c','lV1p','lOFA','lFFA','lATC','lmPFC','lMTG','lAG','lPCC','rFEF','rSPL','rMT','rV1','rV1c','rV1p','rOFA','rFFA','rATC','rMTG','rmPFC','rAG','rPCC'};

% I made the diagonals 0 to threshold so that matrices are not skewed
adj_YC_plot = [1; corr_means.(char(ROIs(1)))(:,1)];
for i=2:length(ROIs)
    adj_YC_plot = [adj_YC_plot [corr_means.(char(ROIs(i)))(1:(i-1),1); 1; corr_means.(char(ROIs(i)))(i:end,1)]];
end

adj_OC_plot = [0; corr_means.(char(ROIs(1)))(:,2)];
for i=2:length(ROIs)
    adj_OC_plot = [adj_OC_plot [corr_means.(char(ROIs(i)))(1:(i-1),2); 1; corr_means.(char(ROIs(i)))(i:end,2)]];
end

adj_svPPA_plot = [0; corr_means.(char(ROIs(1)))(:,3)];
for i=2:length(ROIs)
    adj_svPPA_plot = [adj_svPPA_plot [corr_means.(char(ROIs(i)))(1:(i-1),3); 1; corr_means.(char(ROIs(i)))(i:end,3)]];
end



% Remove unwanted rows and columns
adj_YC_plot(:,6) = [];
adj_YC_plot(:,4) = [];
adj_YC_plot(6,:) = [];
adj_YC_plot(4,:) = [];
adj_YC_plot(:,17) = [];
adj_YC_plot(:,15) = [];
adj_YC_plot(17,:) = [];
adj_YC_plot(15,:) = [];

adj_OC_plot(:,6) = [];
adj_OC_plot(:,4) = [];
adj_OC_plot(6,:) = [];
adj_OC_plot(4,:) = [];
adj_OC_plot(:,17) = [];
adj_OC_plot(:,15) = [];
adj_OC_plot(17,:) = [];
adj_OC_plot(15,:) = [];

adj_svPPA_plot(:,6) = [];
adj_svPPA_plot(:,4) = [];
adj_svPPA_plot(6,:) = [];
adj_svPPA_plot(4,:) = [];
adj_svPPA_plot(:,17) = [];
adj_svPPA_plot(:,15) = [];
adj_svPPA_plot(17,:) = [];
adj_svPPA_plot(15,:) = [];


ROIs = {'lFEF','lSPL','lMT','lV1','lV1c','lV1p','lOFA','lFFA','lATC','lmPFC','lMTG','lAG','lPCC','rFEF','rSPL','rMT','rV1','rV1c','rV1p','rOFA','rFFA','rATC','rMTG','rmPFC','rAG','rPCC'};
ROIs(6) = [];
ROIs(4) = [];
ROIs(17) = [];
ROIs(15) = [];

figure('Renderer', 'painters', 'Position', [10 10 900 600])
imagesc(adj_svPPA_plot);
colorbar
colormap(redblue)
set(gca, 'XTick', 1:22, 'XTickLabel', ROIs, 'YTick', 1:22, 'YTickLabel', ROIs)
set(gca, 'FontSize', 18)
caxis([-1 1])
xtickangle(45)
title('svPPA','fontsize',24);
