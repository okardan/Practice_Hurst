% Reference:
% Kardan, O., et al. (2023) Improvements in task performance after practice
% are associated with scale-free dynamics of brain activity.
% Network Neuroscience, 1-63. https://doi.org/10.1162/netn_a_00319

% Data used in this script are here: https://osf.io/zsxfj

% This script runs PLS analysis on HCP n-back parcel-wise H data (Study 2)
% Results are used in Figure 4 of the paper and corresponding supplementary
% analyses using Craddock 392 parcellation, task-timings-regressed, or WLMF

clear all
addpath(genpath('~\NIFTI_tools')); % Copyright (c) 2009, Jimmy Shen see license.txt inside NIFTI_tools
shen268 = load_untouch_nii('~\NUBE_light\shen_1mm_268_parcellation.nii'); %Shen et al.,(2013) parcellation
cc392 = load_untouch_nii('~\NUBE_light\CC400.nii'); %Craddock et al.,(2012) parcellation
hcpruns = readtable('~/hcp_run_FDfilter.csv'); %data from HCP study ( https://db.humanconnectome.org) and obtained from Kardan et al., (2022)
run2 = hcpruns(string(hcpruns.Run) == 'LR' & hcpruns.Ns == 2,:);
run1 = hcpruns(string(hcpruns.Run) == 'RL' & hcpruns.Ns == 2,:);
hcpdat = innerjoin(run1,run2,'Key','subs');
%%
figure
hold on
for i=1:599
    plot(linspace(0,1,100),linspace(hcpdat.Acc_run1(i),hcpdat.Acc_run2(i)),'color',[.8,.7,.9]);
end
plot(linspace(0,1,100),linspace(mean(hcpdat.Acc_run1),mean(hcpdat.Acc_run2)),'color',[.1,.1,.1]);
hold off

[h,p,ci,stats] = ttest(hcpdat.Acc_run1,hcpdat.Acc_run2)
%%
addpath(genpath('~\PlS')); % download PLS containing plscmd and plsgui folders from https://github.com/McIntosh-Lab/PLS
yy = hcpdat.Acc_run2 - hcpdat.Acc_run1; XX = [ones(599,1) hcpdat.Acc_run1];
[b,bint,classlab,rint,stats]  = regress(yy,XX);
incsubs =1:599;
rng('default');
groups{1} = 1:599;
nRuns = 2;Runs=[1,2];
% switch between Shen parcels or Craddock 392 (CC400) parcels
nParcels = 268; badParcels =[]; goodps = setdiff(1:268,badParcels); shencradd = 'Shen';
% nParcels = 392; badParcels=[]; goodps = setdiff(1:392,badParcels); shencradd = 'Cradd';

inds = find(tril(ones(length(goodps)),-1)==1);
datamat_lst = cell(length(groups),1);
% ####*****************####
delaprimes = classlab; % or raw yy %adjusted change in perf or non-adjusted for supplementary

for g = 1:numel(groups)
    for c = 1:nRuns
        for s = 1:numel(groups{g})
            
            
            
            if c==1
                %                 load(['~\cc400_Hs_studies1-3\hcp_nbk\',num2str(hcpdat.subs(s)),'_RL_cc4_H.mat']);
                %                 vec = voxel_H(:,1)'; %for supplementary Craddock analysis
                
                %                 load(['~\task_regged_Hs_shen\hcp_regged_',num2str(hcpdat.subs(s)),'_RL_shen_H.mat']);
                %                  vec = parcel_H(:,1)';  %for supplementary task_regressed analysis
                
                load(['~\hcp_shenH\',num2str(hcpdat.subs(s)),'_RL_shen_H.mat']);
                %                  vec = shen_H(:,2)'; % c1 for supplementary WLMF analysis
                vec = shen_H(:,1)'; %
            end
            if c==2
                %                 load(['~\cc400_Hs_studies1-3\hcp_nbk\',num2str(hcpdat.subs(s)),'_LR_cc4_H.mat']);
                %                   vec = voxel_H(:,1)'; %for supplementary Craddock analysis
                
                %                 load(['~\task_regged_Hs_shen\hcp_regged_',num2str(hcpdat.subs(s)),'_LR_shen_H.mat']);
                %                 vec = parcel_H(:,1)'; %for supplementary task_regressed analysis
                
                load(['~\hcp_shenH\',num2str(hcpdat.subs(s)),'_LR_shen_H.mat']);
                %                 vec = shen_H(:,2)'; %c1 for supplementary WLMF analysis
                vec = shen_H(:,1)'; %
            end
            
            datamat_lst{g} = [datamat_lst{g}; vec];
            s
        end
    end
end


num_subj = [length(groups{1})];
num_cond = nRuns;

option.method = 3; %behavioral PLS
option.num_boot = 5000;
option.num_perm = 1000;
option.meancentering_type=[2]; % 2 is for only grand mean removed
option.stacked_behavdata = [delaprimes; delaprimes]; %
result = pls_analysis(datamat_lst, num_subj, num_cond, option);

%put in the datamat_lst into results.  No longer automatic in new PLS
%version

result.datamat_lst = datamat_lst;
ps = result.perm_result.sprob
brainPs = result.boot_result.compare_u(result.boot_result.compare_u(:,1)>3)
goodps(find(result.boot_result.compare_u(:,1)>3))
brainNs = result.boot_result.compare_u(result.boot_result.compare_u(:,1)<-3)
goodps(find(result.boot_result.compare_u(:,1)<-3))
figure;
for k=1:nRuns
    scatter(k,result.boot_result.ulcorr_adj(k,1));hold on
    scatter(k,result.boot_result.orig_corr(k,1));
    scatter(k,result.boot_result.llcorr_adj(k,1));
end
hold off

%save('2BKadj_Shen_H_behpls_result.mat','result');

%% beh PLS figure making
%load('2BKadj_Shen_H_behpls_result.mat')
for lv=1:1
    %%
    figure;
    %flipped sign on both sides of LV for ease of interpretation
    line(ones(300,1),linspace(0,-result.boot_result.orig_corr(1,1),300),'Color','r','LineWidth',13);
    hold on
    for k=1:2
        
        col = 'k';
        xpos = k ;
        
        line(xpos*ones(300,1),linspace(0,-result.boot_result.orig_corr(k,lv),300),'Color',col,'LineWidth',13);
        scatter(xpos,-result.boot_result.orig_corr(k,lv),'+k');
        line(xpos*ones(300,1),linspace(-result.boot_result.ulcorr_adj(k,lv),-result.boot_result.llcorr_adj(k,lv),300));
    end
    
    sigfeats = find(abs(result.boot_result.compare_u(:,lv))>3); % flipped the sign above in result.boot so flip in brain fig too
    sigfeatsP = find(result.boot_result.compare_u(:,lv)>3);
    sigfeatsN = find(result.boot_result.compare_u(:,lv)<-3);
    z1percent = length(sigfeats==1);
    result.boot_result.compare_u(sigfeats,lv)
    cbcovs = result.s.^2./sum(result.s.^2);
    title(['LV',num2str(lv),'  \sigma_{XY} = ',num2str(cbcovs(lv)),'  p = ',num2str(ps(lv))]);
    ylim([-1,1.2]); xlim([.5,(nRuns +.5)]);
    ylabel('Correlation with adjusted \Delta Accuracy');
    set(gca,'XTick',1:2,'Xticklabel',{'1st NBK run', '2nd NBK run'},...
        'XtickLabelRotation',45,'FontSize',14);
    hold off
    %%
    sigids = goodps(sigfeats);
    sigidsP = goodps(sigfeatsP);
    sigidsN = goodps(sigfeatsN);
    gen_parc_shen = zeros(61,72,56);
    gen_parc_gord = zeros(61,72,56);
    gen_parc_cc4 = zeros(61,72,56);
    for vx=3:61-3
        for vy=3:72-3
            for vz=3:56-3
                xmm = (vx-30.5)*3.25; vxshen = round(91-xmm);  vxcc4 = round(32-xmm/3);
                ymm = (vy-41.6)*3.25; vyshen = round(ymm+127);  vycc4 = round(ymm/3 +44);
                zmm = (vz-23.3)*3.5; vzshen = round(zmm+73);  vzcc4 = round(zmm/3 + 25);
                gen_parc_shen(vx,vy,vz) = shen268.img(vxshen,vyshen,vzshen);
                gen_parc_cc4(vx,vy,vz) = cc392.img(vxcc4,vycc4,vzcc4);
            end
        end
    end
    gen_parc_shenSigs = gen_parc_shen;
    gen_parc_shenSigs(ismember(gen_parc_shen,sigidsP))=-3;%flipped sign on both sides of LV for ease of interpretation
    gen_parc_shenSigs(ismember(gen_parc_shen,sigidsN))=+3;
    gen_parc_shenSigs(~ismember(gen_parc_shen,sigids) & gen_parc_shen~=0)=0;
    nii = make_nii(gen_parc_shenSigs, [3.25,3.25,3.5], [30.5,41.6,23.3]);
    nii.hdr.dime.datatype =64; % used to make Figure 4 in the paper
    save_nii(nii,['~\behPLS_figures\delta2BKadj_sigHs_',shencradd,'_LV',num2str(lv),'.nii']);
    
    % for the craddock version
    %     gen_parc_cc4Sigs = gen_parc_cc4;
    %     gen_parc_cc4Sigs(ismember(gen_parc_cc4,sigidsP))=+3;%
    %     gen_parc_cc4Sigs(ismember(gen_parc_cc4,sigidsN))=-3;
    %     gen_parc_cc4Sigs(~ismember(gen_parc_cc4,sigids) & gen_parc_cc4~=0)=0;
    %
    %
    %     nii = make_nii(gen_parc_cc4Sigs, [3.25,3.25,3.5], [30.5,41.6,23.3]);
    %     nii.hdr.dime.datatype =64;
    %     save_nii(nii,['~\behPLS_figures\delta2BKadj_sigHs_',shencradd,'_LV',num2str(lv),'.nii']);
    
end