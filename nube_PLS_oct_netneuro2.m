% Reference:
% Kardan, O., et al. (2023) Improvements in task performance after practice
% are associated with scale-free dynamics of brain activity.
% Network Neuroscience, 1-63. https://doi.org/10.1162/netn_a_00319

% Data used in this script are here: https://osf.io/zsxfj

% This script runs PLS analysis on dual n-back data parcel-wise H (Study 1)
% Results are used in Figure 3 of the paper and corresponding supplementary
% analyses using Craddock 392 parcellation, task-timings-regressed, or WLMF

clear all

addpath(genpath('~\NIFTI_tools')); % Copyright (c) 2009, Jimmy Shen see license.txt inside NIFTI_tools
shen268 = load_untouch_nii('~\NUBE_light\shen_1mm_268_parcellation.nii'); %Shen et al.,(2013) parcellation
cc392 = load_untouch_nii('~\NUBE_light\CC400.nii'); %Craddock et al.,(2012) parcellation

addpath(genpath('~\PlS')); % download PLS containing plscmd and plsgui folders from https://github.com/McIntosh-Lab/PLS
load('~\Shen_Parcel_H_data_studies1&3\NubeOct2023.mat')
Both_diff1 =[]; baseA=[];
for i=1:56
    baseA = [baseA; NubeOct2023{i, 1}.DNB1.Aprm];
    Both_diff1 = [Both_diff1; NubeOct2023{i, 1}.DNB2.Aprm - NubeOct2023{i, 1}.DNB1.Aprm];
end

%incsubs =[1:41,43:56]; shencradd = 'Cradd';% for cc-400 parcellation one sub was removed as they had 13 NaN parcels
incsubs =[1:56]; shencradd = 'Shen';

yy = Both_diff1(incsubs); XX = [ones(length(incsubs),1) baseA(incsubs)];
[b,bint,classlab,rint,stats]  = regress(yy,XX);


rng('default');

groups{1} = incsubs;
nRuns = 3;Runs=[1,2,3];

% switch between Shen parcels or Craddock 392 (CC400) parcels
nParcels = 268; badParcels =[4,131,137,189,239,252]; goodps = setdiff(1:268,badParcels); shencradd = 'Shen'; %parcels with NaN values in Study 1 and 3 are excluded
% nParcels = 392; badParcels=[48]; goodps = setdiff(1:392,badParcels); shencradd = 'Cradd';

%%
datamat_lst = cell(length(groups),1); 
%#####********
delaprimes = classlab; % or yy %adjusted change in perf or non-adj for supp

for g = 1:numel(groups)
    for c = 1:nRuns
        for s = 1:numel(groups{g})
         ss = incsubs(s);   
            
%             hmat1 = NubeOct2023{ss, 1}.DNB1.craddH(goodps,1)';  %for supplementary Craddock analysis
%             hmat2 = NubeOct2023{ss, 1}.vid.craddH(goodps,1)';
%             hmat3 = NubeOct2023{ss, 1}.DNB2.craddH(goodps,1)';
%             

%             hmat1 = NubeOct2023{ss, 1}.DNB1.tskreggedshenH(goodps,1)';  %for supplementary task_regressed analysis
%             hmat2 = NubeOct2023{ss, 1}.vid.shH(goodps,1)';
%             hmat3 = NubeOct2023{ss, 1}.DNB2.tskreggedshenH(goodps,1)';
%             
            
%             hmat1 = NubeOct2023{ss, 1}.DNB1.shH(goodps,2)'; % for WLMF first cumulant
%             hmat2 = NubeOct2023{ss, 1}.vid.shH(goodps,2)';
%             hmat3 = NubeOct2023{ss, 1}.DNB2.shH(goodps,2)';
%             

%             hmat1 = NubeOct2023{ss, 1}.DNB1.shH(goodps,3)'; %2nd cumulant
%             hmat2 = NubeOct2023{ss, 1}.vid.shH(goodps,3)';
%             hmat3 = NubeOct2023{ss, 1}.DNB2.shH(goodps,3)';
%             

            hmat1 = NubeOct2023{ss, 1}.DNB1.shH(goodps,1)'; % main analysis
            hmat2 = NubeOct2023{ss, 1}.vid.shH(goodps,1)';
            hmat3 = NubeOct2023{ss, 1}.DNB2.shH(goodps,1)';
%           
            
            if c ==1
                vec = hmat1;
            end
            if c ==2
                vec = hmat2;
            end
            if c ==3
                vec = hmat3;
            end
            datamat_lst{g} = [datamat_lst{g}; vec];
            s
        end
    end
end


% num_subj = [length(groups{1})  length(groups{2})];
num_subj = [length(groups{1})];
num_cond = nRuns;

option.method = 3; %behavioral PLS 3
option.num_boot = 5000;
option.num_perm = 1000;
option.meancentering_type=[2];
option.cormode = 0; %	0. Pearson correlation 2. covaraince 4. cosine angle 6. dot product
option.stacked_behavdata = [delaprimes; delaprimes; delaprimes;]
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

%%
%save('NUBEadj_Shen_H_behpls_result.mat','result');

%% beh PLS figure making
%load('NUBEadj_Shen_H_behpls_result.mat')
for lv=1:1
    %%
    figure;
    %flipped sign on both sides of LV for ease of interpretation
    line(ones(300,1),linspace(0,-result.boot_result.orig_corr(1,1),300),'Color','r','LineWidth',13);
     hold on
for k=1:3
   
col = 'k';
    xpos = k ;
    line(xpos*ones(300,1),linspace(0,-result.boot_result.orig_corr(k,lv),300),'Color',col,'LineWidth',13);
    scatter(xpos,-result.boot_result.orig_corr(k,lv),'+k');
    line(xpos*ones(300,1),linspace(-result.boot_result.ulcorr_adj(k,lv),-result.boot_result.llcorr_adj(k,lv),300));
end
    
    sigfeats = find(abs(result.boot_result.compare_u(:,lv))>2); % flipped the sign above in result.boot so flip in brain fig too
    sigfeatsP = find(result.boot_result.compare_u(:,lv)>2);
    sigfeatsN = find(result.boot_result.compare_u(:,lv)<-2);
    z1percent = length(sigfeats==1);
    result.boot_result.compare_u(sigfeats,lv)
    cbcovs = result.s.^2./sum(result.s.^2);
    title(['LV',num2str(lv),'  \sigma_{XY} = ',num2str(cbcovs(lv)),'  p = ',num2str(ps(lv))]);
    ylim([-1,1.2]); xlim([.5,(nRuns +.5)]);
    ylabel('Correlation with adjusted (\Delta A'')');
    set(gca,'XTick',1:nRuns,...
        'Xticklabel',{'1st DNB run','Video','2nd DNB run'},...
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
    gen_parc_shenSigs(ismember(gen_parc_shen,sigidsP))=-2;%flipped sign on both sides of LV for ease of interpretation
    gen_parc_shenSigs(ismember(gen_parc_shen,sigidsN))=+2;
    gen_parc_shenSigs(~ismember(gen_parc_shen,sigids) & gen_parc_shen~=0)=0;
    
    nii = make_nii(gen_parc_shenSigs, [3.25,3.25,3.5], [30.5,41.6,23.3]);
    nii.hdr.dime.datatype =64;
    save_nii(nii,['~\behPLS_figures\deltaDNBregBase_sigzHs_',shencradd,'_LV',num2str(lv),'.nii']);
    
    % for the craddock version
%     gen_parc_cc4Sigs = gen_parc_cc4;
%     gen_parc_cc4Sigs(ismember(gen_parc_cc4,sigidsP))=+2;%
%     gen_parc_cc4Sigs(ismember(gen_parc_cc4,sigidsN))=-2;
%     gen_parc_cc4Sigs(~ismember(gen_parc_cc4,sigids) & gen_parc_cc4~=0)=0;
%     
%     nii = make_nii(gen_parc_cc4Sigs, [3.25,3.25,3.5], [30.5,41.6,23.3]);
%     nii.hdr.dime.datatype =64;
%     save_nii(nii,['~\behPLS_figures\deltaDNBregBase_sigzHs_',shencradd,'_LV',num2str(lv),'.nii']);

end