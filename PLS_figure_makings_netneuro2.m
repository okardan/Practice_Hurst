% Reference:
% Kardan, O., et al. (2023) Improvements in task performance after practice
% are associated with scale-free dynamics of brain activity.
% Network Neuroscience, 1-63. https://doi.org/10.1162/netn_a_00319

% Data used in this script are here: https://osf.io/zsxfj and also produced
% by the other scripts that perform PLs analysis

% This script uses the PLS analysis results on each datsaset to make
% the behavior-brain plots that are tab iv in Figures 3-5 and thier
% supplementary variants (Craddock 392 parcellation,
% task-timings-regressed, or WLMF)


%% dual n-back (study 1, fig 3 and supp variants)
clear all

addpath(genpath('~\NIFTI_tools')); % Copyright (c) 2009, Jimmy Shen see license.txt inside NIFTI_tools
shen268 = load_untouch_nii('~\NUBE_light\shen_1mm_268_parcellation.nii'); %Shen et al.,(2013) parcellation
cc392 = load_untouch_nii('~\NUBE_light\CC400.nii'); %Craddock et al.,(2012) parcellation
load('~\Shen_Parcel_H_data_studies1&3\NubeOct2023.mat')
Both_diff1 =[]; baseA=[];
for i=1:56
    baseA = [baseA; NubeOct2023{i, 1}.DNB1.Aprm];
    Both_diff1 = [Both_diff1; NubeOct2023{i, 1}.DNB2.Aprm - NubeOct2023{i, 1}.DNB1.Aprm];
end
% incsubs =[1:41,43:56]; shencradd = 'Cradd';% for cc-400 parcellation one sub was removed as they had 13 NaN parcels
incsubs =[1:56]; shencradd = 'Shen';
yy = Both_diff1(incsubs); XX = [ones(length(incsubs),1) baseA(incsubs)];
[b,bint,classlab,rint,stats]  = regress(yy,XX);
rng('default');

groups{1} = incsubs;
nRuns = 3;Runs=[1,2,3];

nParcels = 268; badParcels =[4,131,137,189,239,252]; goodps = setdiff(1:268,badParcels); shencradd = 'Shen';
% nParcels = 392; badParcels=[48]; goodps = setdiff(1:392,badParcels); shencradd = 'Cradd';

% PLS results from the the other scripts; choose one to be used in the figure (top line is main analysis, others are supp) 
 load('~\behPLS_figures\NUBEadj_Shen_H_behpls_result.mat')
%  load('~\behPLS_figures\NUBEadj_Shen_H_behplscc4_result.mat')
%   load('~\behPLS_figures\NUBEadj_Shen_H_behplstaskregged_result.mat')
%   load('~\behPLS_figures\NUBEadj_Shen_H_behplsWLMF_result.mat')
%   load('~\behPLS_figures\NUBEnonadj_Shen_H_behpls_result.mat')
  
[length(goodps) length(result.u)]

lv=1;
    sigfeats = find(abs(result.boot_result.compare_u(:,lv))>3); % 99% CI
    sigfeatsP = find(result.boot_result.compare_u(:,lv)>3);
    sigfeatsN = find(result.boot_result.compare_u(:,lv)<-3);
    
    sigids = goodps(sigfeats);
    sigidsP = goodps(sigfeatsP);
    sigidsN = goodps(sigfeatsN);
    gen_parc_shen = zeros(61,72,56);
    gen_parc_gord = zeros(61,72,56);
    gen_parc_cc4 = zeros(61,72,56);
    for vx=3:61-3
        for vy=3:72-3
            for vz=3:56-3
                xmm = (vx-30.5)*3.25; vxshen = round(91-xmm); vxgord = round(xmm+92); vxcc4 = round(32-xmm/3);
                ymm = (vy-41.6)*3.25; vyshen = round(ymm+127); vygord = round(ymm+127); vycc4 = round(ymm/3 +44);
                zmm = (vz-23.3)*3.5; vzshen = round(zmm+73); vzgord = round(zmm+73); vzcc4 = round(zmm/3 + 25);
                gen_parc_shen(vx,vy,vz) = shen268.img(vxshen,vyshen,vzshen);
                gen_parc_gord(vx,vy,vz) = gordon333.img(vxgord,vygord,vzgord);
                gen_parc_cc4(vx,vy,vz) = cc392.img(vxcc4,vycc4,vzcc4);
            end
        end
    end
    gen_parc_shenSigs = gen_parc_shen;
    gen_parc_shenSigs(ismember(gen_parc_shen,sigidsP))=-3;%flipped sign on both sides of LV for ease of interpretation
    gen_parc_shenSigs(ismember(gen_parc_shen,sigidsN))=+3;
    gen_parc_shenSigs(~ismember(gen_parc_shen,sigids) & gen_parc_shen~=0)=0;
    nii = make_nii(gen_parc_shenSigs, [3.25,3.25,3.5], [30.5,41.6,23.3]);
    nii.hdr.dime.datatype =64;
    save_nii(nii,['~\behPLS_figures\deltaDNBadj_sigHs_',shencradd,'_LV',num2str(lv),'.nii']);
    
%     gen_parc_cc4Sigs = gen_parc_cc4;
%     gen_parc_cc4Sigs(ismember(gen_parc_cc4,sigidsP))=+3;%
%     gen_parc_cc4Sigs(ismember(gen_parc_cc4,sigidsN))=-3;
%     gen_parc_cc4Sigs(~ismember(gen_parc_cc4,sigids) & gen_parc_cc4~=0)=0;
%   
%     nii = make_nii(gen_parc_cc4Sigs, [3.25,3.25,3.5], [30.5,41.6,23.3]);
%     nii.hdr.dime.datatype =64;
%     save_nii(nii,['~\behPLS_figures\deltaDNBregBase_sigHs_',shencradd,'_LV',num2str(lv),'.nii']);


% NUBE figure beh pls
 [r,p] = corr(result.usc(:,1), result.vsc(:,1)); r^2
sg = +1; %flip sign if needed
 ps = result.perm_result.sprob
 brainPs = result.boot_result.compare_u(result.boot_result.compare_u(:,1)>3)
goodps(find(result.boot_result.compare_u(:,1)>3))
brainNs = result.boot_result.compare_u(result.boot_result.compare_u(:,1)<-3)
goodps(find(result.boot_result.compare_u(:,1)<-3))
% 
for lv=1:1
    figure('InnerPosition',[209,601,795,457],'OuterPosition',[201,593,811,550]) ;
    err = [result.boot_result.llcorr_adj(:,lv)',result.boot_result.ulcorr_adj(:,lv)'];
    bar([1,2,3],sg*result.lvcorrs(:,lv)',.4,'FaceColor', [0.96,0.78,0.35], 'EdgeColor','none');hold on
    if sg == -1
    errorbar([1,2,3],sg*result.lvcorrs(:,lv)',...
            sg*(result.boot_result.ulcorr_adj(:,lv)' - result.lvcorrs(:,lv)'),...
            sg*(result.lvcorrs(:,lv)'-result.boot_result.llcorr_adj(:,lv)'),'+','color',[0.07,0.62,1.00],'CapSize',0); 
    end
    if sg== +1
    errorbar([1,2,3],sg*result.lvcorrs(:,lv)',...
             sg*(result.lvcorrs(:,lv)'-result.boot_result.llcorr_adj(:,lv)'),...
           sg*(result.boot_result.ulcorr_adj(:,lv)' - result.lvcorrs(:,lv)'),'+','color',[0.07,0.62,1.00],'CapSize',0);         
    end   
    sigfeats = find(abs(result.boot_result.compare_u(:,lv))>3); % 99% CI
    sigfeatsP = find(result.boot_result.compare_u(:,lv)>3);
    sigfeatsN = find(result.boot_result.compare_u(:,lv)<-3);
    z1percent = length(sigfeats==1)
    result.boot_result.compare_u(sigfeats,lv)
    cbcovs = result.s.^2./sum(result.s.^2);
%     title(['LV',num2str(lv),'  \sigma_{XY} = ',num2str(cbcovs(lv)),'  p = ',num2str(round(ps(lv)*1000)/1000)]);
disp(['LV',num2str(lv),'  \sigma_{XY} = ',num2str(cbcovs(lv)),'  p = ',num2str(round(ps(lv)*1000)/1000)])
    ylim([-1,1]); %xlim([.9,2.1]);
    ylabel('Correlation with adjusted \Delta A''');
%     ylabel('Correlation with \Delta A''');
    set(gca,'XTick',1:3,'Xticklabel',{'', '', ''},...
        'XtickLabelRotation',45,'FontSize',14);
    text(.8,-.9,'1st DNB','Rotation',45,'FontSize',14); text(1.8,-.9,'video','Rotation',45,'FontSize',14); text(2.8,-.9,'2nd DNB','Rotation',45,'FontSize',14)
    hold off
end
%% save figures using export_fig: Yair Altman (2023). export_fig (https://github.com/altmany/export_fig/releases/tag/v3.37
% addpath(genpath('~\github_repo_export_fig'))
% export_fig myfigurename.jpg -m 4 -transparent
%% n-back (study 2, fig 4 and supp variants)
clear all

addpath(genpath('~\NIFTI_tools')); % Copyright (c) 2009, Jimmy Shen see license.txt inside NIFTI_tools
shen268 = load_untouch_nii('~\NUBE_light\shen_1mm_268_parcellation.nii'); %Shen et al.,(2013) parcellation
cc392 = load_untouch_nii('~\NUBE_light\CC400.nii'); %Craddock et al.,(2012) parcellation

nParcels = 268; badParcels =[]; goodps = setdiff(1:268,badParcels); shencradd = 'Shen';
% nParcels = 392; badParcels=[]; goodps = setdiff(1:392,badParcels); shencradd = 'Cradd';

% PLS results from the the other scripts; choose one to be used in the figure (top line is main analysis, others are supp)
 load('~\behPLS_figures\2BKadj_Shen_H_behpls_result.mat')
%  load('~\behPLS_figures\2BKadj_Shen_H_behplscc4_result.mat')
%  load('~\behPLS_figures\2BKadj_Shen_H_behplstaskregged_result.mat')
% load('~\behPLS_figures\2BKadj_Shen_H_behplsWLMF_result.mat')
% load('~\behPLS_figures\2BKnonadj_Shen_H_behpls_result.mat')

[length(goodps) length(result.u)]

lv=1;
    sigfeats = find(abs(result.boot_result.compare_u(:,lv))>3); % 99% CI
    sigfeatsP = find(result.boot_result.compare_u(:,lv)>3);
    sigfeatsN = find(result.boot_result.compare_u(:,lv)<-3);

    sigids = goodps(sigfeats);
    sigidsP = goodps(sigfeatsP);
    sigidsN = goodps(sigfeatsN);
    gen_parc_shen = zeros(61,72,56);
    gen_parc_gord = zeros(61,72,56);
    gen_parc_cc4 = zeros(61,72,56);
    for vx=3:61-3
        for vy=3:72-3
            for vz=3:56-3
                xmm = (vx-30.5)*3.25; vxshen = round(91-xmm); vxgord = round(xmm+92); vxcc4 = round(32-xmm/3);
                ymm = (vy-41.6)*3.25; vyshen = round(ymm+127); vygord = round(ymm+127); vycc4 = round(ymm/3 +44);
                zmm = (vz-23.3)*3.5; vzshen = round(zmm+73); vzgord = round(zmm+73); vzcc4 = round(zmm/3 + 25);
                gen_parc_shen(vx,vy,vz) = shen268.img(vxshen,vyshen,vzshen);
                gen_parc_gord(vx,vy,vz) = gordon333.img(vxgord,vygord,vzgord);
                gen_parc_cc4(vx,vy,vz) = cc392.img(vxcc4,vycc4,vzcc4);
            end
        end
    end
    gen_parc_shenSigs = gen_parc_shen;
    sg = -1; %flip sign if needed
    gen_parc_shenSigs(ismember(gen_parc_shen,sigidsP))=sg*3;%flipped sign on both sides of LV for ease of interpretation
    gen_parc_shenSigs(ismember(gen_parc_shen,sigidsN))=sg*(-3);
    gen_parc_shenSigs(~ismember(gen_parc_shen,sigids) & gen_parc_shen~=0)=0;
 
    nii = make_nii(gen_parc_shenSigs, [3.25,3.25,3.5], [30.5,41.6,23.3]);
    nii.hdr.dime.datatype =64;
    save_nii(nii,['~\behPLS_figures\delta2BKadj_sigHs_',shencradd,'_LV',num2str(lv),'.nii']);
  

%     gen_parc_cc4Sigs = gen_parc_cc4;
%     gen_parc_cc4Sigs(ismember(gen_parc_cc4,sigidsP))=+3;%
%     gen_parc_cc4Sigs(ismember(gen_parc_cc4,sigidsN))=-3;
%     gen_parc_cc4Sigs(~ismember(gen_parc_cc4,sigids) & gen_parc_cc4~=0)=0;
    
  
%     nii = make_nii(gen_parc_cc4Sigs, [3.25,3.25,3.5], [30.5,41.6,23.3]);
%     nii.hdr.dime.datatype =64;
%     save_nii(nii,['~\behPLS_figures\delta2BKadj_sigHs_',shencradd,'_LV',num2str(lv),'.nii']);

   % % behavioral loadings
    [r,p] = corr(result.usc(:,1), result.vsc(:,1)); r^2

 ps = result.perm_result.sprob
 brainPs = result.boot_result.compare_u(result.boot_result.compare_u(:,1)>3)
goodps(find(result.boot_result.compare_u(:,1)>3))
brainNs = result.boot_result.compare_u(result.boot_result.compare_u(:,1)<-3)
goodps(find(result.boot_result.compare_u(:,1)<-3))
% 
for lv=1:1
    figure('InnerPosition',[209,601,795,457],'OuterPosition',[201,593,811,550]) ;
    err = [result.boot_result.llcorr_adj(:,lv)',result.boot_result.ulcorr_adj(:,lv)'];
    bar([1,2],sg*result.lvcorrs(:,lv)',.4,'FaceColor', [0.96,0.78,0.35], 'EdgeColor','none');hold on
    if sg == -1
    errorbar([1,2],sg*result.lvcorrs(:,lv)',...
            sg*(result.boot_result.ulcorr_adj(:,lv)' - result.lvcorrs(:,lv)'),...
            sg*(result.lvcorrs(:,lv)'-result.boot_result.llcorr_adj(:,lv)'),'+','color',[0.07,0.62,1.00],'CapSize',0); 
    end
    if sg== +1
    errorbar([1,2],sg*result.lvcorrs(:,lv)',...
             sg*(result.lvcorrs(:,lv)'-result.boot_result.llcorr_adj(:,lv)'),...
           sg*(result.boot_result.ulcorr_adj(:,lv)' - result.lvcorrs(:,lv)'),'+','color',[0.07,0.62,1.00],'CapSize',0);         
    end   
    sigfeats = find(abs(result.boot_result.compare_u(:,lv))>3); % 99% CI
    sigfeatsP = find(result.boot_result.compare_u(:,lv)>3);
    sigfeatsN = find(result.boot_result.compare_u(:,lv)<-3);
    z1percent = length(sigfeats==1)
    result.boot_result.compare_u(sigfeats,lv)
    cbcovs = result.s.^2./sum(result.s.^2);
%     title(['LV',num2str(lv),'  \sigma_{XY} = ',num2str(cbcovs(lv)),'  p = ',num2str(round(ps(lv)*1000)/1000)]);
disp(['LV',num2str(lv),'  \sigma_{XY} = ',num2str(cbcovs(lv)),'  p = ',num2str(round(ps(lv)*1000)/1000)])
    ylim([-1,1]); %xlim([.9,2.1]);
    ylabel('Correlation with adjusted \Delta Accuracy');
%     ylabel('Correlation with \Delta Accuracy');
    set(gca,'XTick',1:3,'Xticklabel',{'', '', ''},...
        'XtickLabelRotation',45,'FontSize',14);
    text(.8,-.9,'1st NBK','Rotation',45,'FontSize',14); text(1.8,-.9,'2nd NBK','Rotation',45,'FontSize',14); 
       hold off
end
%% save figures using export_fig: Yair Altman (2023). export_fig (https://github.com/altmany/export_fig/releases/tag/v3.37
% addpath(genpath('~\github_repo_export_fig'))
% export_fig myfigurename.jpg -m 4 -transparent
%% choose-and-solve task (study 3, fig 5 and supp variants)
clear all

addpath(genpath('E:\Omid\MatlabGraphics\NIFTI_tools'));
shen268 = load_untouch_nii('E:\Omid\MatlabGraphics\NUBE_light\shen_1mm_268_parcellation.nii');
gordon333 = load_untouch_nii('E:\Omid\MatlabGraphics\NUBE_light\gordon333MNI.nii');
cc392 = load_untouch_nii('E:\Omid\MatlabGraphics\NUBE_light\CC400.nii');

nParcels = 268; badParcels =[4,131,137,189,239,252]; goodps = setdiff(1:268,badParcels); shencradd = 'Shen';
% nParcels = 392; badParcels=[48]; goodps = setdiff(1:392,badParcels); shencradd = 'Cradd';

% PLS results from the the other scripts; choose one to be used in the figure (top line is main analysis, others are supp)
 load('~\behPLS_figures\CASTadj_Shen_H_behpls_result.mat')
%  load('~\behPLS_figures\CASTadj_Shen_H_behplscc4_result.mat')
%  load('~\behPLS_figures\CASTadj_Shen_H_behplstaskregged_result.mat')
%  load('~\behPLS_figures\CASTadj_Shen_H_behplsWLMF_result.mat')
%  load('~\behPLS_figures\CASTnonadj_Shen_H_behpls_result.mat')
%  load('~\behPLS_figures\CASTmeancenteredwithin_Shen_H_behpls_result.mat')
% load('~\behPLS_figures\CASTmeancenteredwithin2_Shen_H_behpls_result.mat')

[length(goodps) length(result.u)]

lv=1;
    sigfeats = find(abs(result.boot_result.compare_u(:,lv))>3); % 99% CI
    sigfeatsP = find(result.boot_result.compare_u(:,lv)>3);
    sigfeatsN = find(result.boot_result.compare_u(:,lv)<-3);

    sigids = goodps(sigfeats);
    sigidsP = goodps(sigfeatsP);
    sigidsN = goodps(sigfeatsN);
    gen_parc_shen = zeros(61,72,56);
    gen_parc_gord = zeros(61,72,56);
    gen_parc_cc4 = zeros(61,72,56);
    for vx=3:61-3
        for vy=3:72-3
            for vz=3:56-3
                xmm = (vx-30.5)*3.25; vxshen = round(91-xmm); vxgord = round(xmm+92); vxcc4 = round(32-xmm/3);
                ymm = (vy-41.6)*3.25; vyshen = round(ymm+127); vygord = round(ymm+127); vycc4 = round(ymm/3 +44);
                zmm = (vz-23.3)*3.5; vzshen = round(zmm+73); vzgord = round(zmm+73); vzcc4 = round(zmm/3 + 25);
                gen_parc_shen(vx,vy,vz) = shen268.img(vxshen,vyshen,vzshen);
                gen_parc_gord(vx,vy,vz) = gordon333.img(vxgord,vygord,vzgord);
                gen_parc_cc4(vx,vy,vz) = cc392.img(vxcc4,vycc4,vzcc4);
            end
        end
    end
    gen_parc_shenSigs = gen_parc_shen;
    sg = +1; %flip sign if needed
    gen_parc_shenSigs(ismember(gen_parc_shen,sigidsP))=sg*3;%flipped sign on both sides of LV for ease of interpretation
    gen_parc_shenSigs(ismember(gen_parc_shen,sigidsN))=sg*(-3);
    gen_parc_shenSigs(~ismember(gen_parc_shen,sigids) & gen_parc_shen~=0)=0;
 
    nii = make_nii(gen_parc_shenSigs, [3.25,3.25,3.5], [30.5,41.6,23.3]);
    nii.hdr.dime.datatype =64;
    save_nii(nii,['~\behPLS_figures\deltaCASTadj_sigHs_',shencradd,'_LV',num2str(lv),'.nii']);
   
%     gen_parc_cc4Sigs = gen_parc_cc4;
%     gen_parc_cc4Sigs(ismember(gen_parc_cc4,sigidsP))=+3;%
%     gen_parc_cc4Sigs(ismember(gen_parc_cc4,sigidsN))=-3;
%     gen_parc_cc4Sigs(~ismember(gen_parc_cc4,sigids) & gen_parc_cc4~=0)=0;
    
 
%     nii = make_nii(gen_parc_cc4Sigs, [3.25,3.25,3.5], [30.5,41.6,23.3]);
%     nii.hdr.dime.datatype =64;
%     save_nii(nii,['E:\Omid\nube\October\Network Neurosciecne\Net Neuro Revisions\behPLS_figures\deltaCASTadj_sigHs_',shencradd,'_LV',num2str(lv),'.nii']);

   % % behavioral loadings
    [r,p] = corr(result.usc(:,1), result.vsc(:,1)); r^2

 ps = result.perm_result.sprob
 brainPs = result.boot_result.compare_u(result.boot_result.compare_u(:,1)>3)
goodps(find(result.boot_result.compare_u(:,1)>3))
brainNs = result.boot_result.compare_u(result.boot_result.compare_u(:,1)<-3)
goodps(find(result.boot_result.compare_u(:,1)<-3))
% 
for lv=1:1
    figure('InnerPosition',[209,601,795,457],'OuterPosition',[201,593,811,550]) ;
    err = [result.boot_result.llcorr_adj(:,lv)',result.boot_result.ulcorr_adj(:,lv)'];
    bar([1:6],sg*result.lvcorrs(:,lv)',.4,'FaceColor', [0.96,0.78,0.35], 'EdgeColor','none');hold on
    if sg == -1
    errorbar([1:6],sg*result.lvcorrs(:,lv)',...
            sg*(result.boot_result.ulcorr_adj(:,lv)' - result.lvcorrs(:,lv)'),...
            sg*(result.lvcorrs(:,lv)'-result.boot_result.llcorr_adj(:,lv)'),'+','color',[0.07,0.62,1.00],'CapSize',0); 
    end
    if sg== +1
    errorbar([1:6],sg*result.lvcorrs(:,lv)',...
             sg*(result.lvcorrs(:,lv)'-result.boot_result.llcorr_adj(:,lv)'),...
           sg*(result.boot_result.ulcorr_adj(:,lv)' - result.lvcorrs(:,lv)'),'+','color',[0.07,0.62,1.00],'CapSize',0);         
    end   
    sigfeats = find(abs(result.boot_result.compare_u(:,lv))>3); % 99% CI
    sigfeatsP = find(result.boot_result.compare_u(:,lv)>3);
    sigfeatsN = find(result.boot_result.compare_u(:,lv)<-3);
    z1percent = length(sigfeats==1)
    result.boot_result.compare_u(sigfeats,lv)
    cbcovs = result.s.^2./sum(result.s.^2);
%     title(['LV',num2str(lv),'  \sigma_{XY} = ',num2str(cbcovs(lv)),'  p = ',num2str(round(ps(lv)*1000)/1000)]);
disp(['LV',num2str(lv),'  \sigma_{XY} = ',num2str(cbcovs(lv)),'  p = ',num2str(round(ps(lv)*1000)/1000)])
    ylim([-1,1]); %xlim([.9,2.1]);
    ylabel('Correlation with adjusted \Delta Accuracy');
%   ylabel('Correlation with \Delta Accuracy');
%   ylabel('Correlation with mean-centered Accuracy');
    set(gca,'XTick',1:6,'Xticklabel',{'', '', '','','',''},...
        'XtickLabelRotation',45,'FontSize',14);
    text(.8,-.9,'1st CAST','Rotation',45,'FontSize',14); text(1.8,-.9,'2nd CAST','Rotation',45,'FontSize',14); 
    text(2.8,-.9,'3rd CAST','Rotation',45,'FontSize',14); text(3.8,-.9,'4th CAST','Rotation',45,'FontSize',14); 
    text(4.8,-.9,'5th CAST','Rotation',45,'FontSize',14); text(5.8,-.9,'6th CAST','Rotation',45,'FontSize',14); 
       hold off
end

%% save figures using export_fig: Yair Altman (2023). export_fig (https://github.com/altmany/export_fig/releases/tag/v3.37
% addpath(genpath('~\github_repo_export_fig'))
% export_fig myfigurename.jpg -m 4 -transparent