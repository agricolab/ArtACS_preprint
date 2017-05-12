matlabrc 
cd('C:\PROJECTS\Global_Projects\ARtACS\')%
addpath('.\repo\src')
addpath('.')
printfolder = '.\docs\img\';

blue    = [0 0.4470  0.7410];
red = [0.8500 0.3250 0.0980];
%%
load('.\data\median\SEP_10Hz.mat');
close all
figure
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')

tmp = trials(2,:);
tmp = reshape(tmp(247:end-254),500,[]);
tmp = tmp-mean(mean(tmp));

xaxis = rad2deg(2*pi*[0.002:0.002:1]);
pick  = 1:300;
hold on
h1 = plot(xaxis(pick),tmp(pick,:),'color',[.5 .5 .5],'linewidth',1);
h2 = plot(xaxis(pick),mean(tmp(pick,:),2),'color',red,'linewidth',2);
set(gca,'XTICK',0:15:360,'XLIM',[60 120])
set(gca,'YTICK',-5000:100:5000,'YLIM',[2750 3300],'YTICKLABEL',[])
xlabel('Phase in degree (°)')
ylabel('Amplitude in a.u.')
l1 = legend([h1(1),h2],'Periods','Average','location','northoutside');
set(l1,'position',[0.15 .85 0 .130])
%title('Distortion')
print(gcf,[printfolder,'intro/distortion.png'],'-dpng')

figure
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
tmp = trials(1,:);
tmp = reshape(tmp(379:end-122),500,[]);
tmp = tmp-mean(mean(tmp));

xaxis = rad2deg(2*pi*[0.002:0.002:1]);
hold on
h1 = plot(xaxis(pick),tmp(pick,:),'color',[.5 .5 .5],'linewidth',1);
h2 = plot(xaxis(pick),mean(tmp(pick,:),2),'color',red,'linewidth',2);
set(gca,'XTICK',0:15:360,'XLIM',[60 120])
set(gca,'YTICK',-5000:100:5000,'YLIM',[2750 3300],'YTICKLABEL',[])
xlabel('Phase in degree (°)')
l1 = legend([h1(1),h2],'Periods','Average','location','northoutside');
set(l1,'position',[0.85 .85 0 .130])
%title('Saturation')
print(gcf,[printfolder,'intro/saturation.png'],'-dpng')
%% ECG PLOTS
NumberPeriods   = 10;
Fs              = 1000;
toi             = 3750:4251;
blperiod        = 3500:3800;
tacsFreq        = 11;
Fs              = 1000;

clean_trials    = utils.getECG('clean.vhdr');
E               = mean(clean_trials,1)-mean(mean(clean_trials,1));
stim_trials     = utils.getECG('11Hz.vhdr');

close all
figure
set(gcf,'Position',[100 100 500 400],'paperpositionmode','auto')
subplot(2,1,1)
hold on
h1 = plot(stim_trials(1,toi),'color',[.8 .8 .8],'linewidth',2);
h2 = plot(mean(clean_trials(:,toi),1),'color',red);
set(gca,'xlim',[1 length(toi)],'xtick',[1:50:4000],'xticklabel',[(-length(toi)/2+1):50:length(toi)/2])
set(gca,'ylim',[-20000 20000],'ytick',[-20000,-10000,0,10000,20000],'yticklabel',[-20,10,0,10,20])
ylabel('mV')
subplot(2,1,2)
hold on
h1 = plot(stim_trials(1,toi),'color',[.8 .8 .8],'linewidth',2);
h2 = plot(E(toi),'r');
set(gca,'xlim',[1 length(toi)],'xtick',[1:50:4000],'xticklabel',[(-length(toi)/2+1):50:length(toi)/2])
set(gca,'ylim',[-20 20],'ytick',[-1000:10:1000])
ylabel('µV')
xlabel('ms')
lh = legend([h1,h2],'Artifacted Recording','Clean Recording');
set(lh,'position',[0.8 .85 .08 .08])
print(gcf,[printfolder,'eva\ecg_raw.png'],'-dpng')

close all
figure
set(gcf,'Position',[100 100 500 200],'paperpositionmode','auto')
hold on
h1 = plot(stim_trials(1,toi),'color',[.8 .8 .8],'linewidth',2);
h2 = plot(mean(clean_trials(:,toi),1),'color',red);
set(gca,'xlim',[1 length(toi)],'xtick',[1:50:4000],'xticklabel',[(-length(toi)/2+1):50:length(toi)/2])
set(gca,'ylim',[-20000 20000],'ytick',[-20000,-10000,0,10000,20000],'yticklabel',[-20,10,0,10,20])
ylabel('mV')
xlabel('ms')
print(gcf,[printfolder,'eva\ecg_raw_1.png'],'-dpng')

figure
set(gcf,'Position',[100 100 500 200],'paperpositionmode','auto')
hold on
h1 = plot(stim_trials(1,toi),'color',[.8 .8 .8],'linewidth',2);
h2 = plot(mean(clean_trials(:,toi),1),'color',red);
set(gca,'xlim',[1 length(toi)],'xtick',[1:50:4000],'xticklabel',[(-length(toi)/2+1):50:length(toi)/2])
set(gca,'ylim',[-20 20],'ytick',[-1000:10:1000])
ylabel('µV')
xlabel('ms')
lh = legend([h1,h2],'Artifacted Recording','Clean Recording');
set(lh,'position',[0.7 .7 0.15 .25],'fontsize',12)
print(gcf,[printfolder,'eva\ecg_raw_2.png'],'-dpng')


R               = [];    
F               = [];

for trl_idx = 1 : size(stim_trials,1)        
    r = stim_trials(trl_idx,:);

    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'causal','uniform','default','default');
    [R(1,trl_idx),F(1,trl_idx,:)]           = utils.regain(f,E,toi,blperiod);
    
    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'causal','linear','default','default');
    [R(2,trl_idx),F(2,trl_idx,:)]           = utils.regain(f,E,toi,blperiod);
    
    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'causal','exp','default','default');
    [R(3,trl_idx),F(3,trl_idx,:)]           = utils.regain(f,E,toi,blperiod);
    
    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'causal','gauss','default','default');
    [R(4,trl_idx),F(4,trl_idx,:)]           = utils.regain(f,E,toi,blperiod);
    
    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'symmetric','automatic','default','default');
    [R(5,trl_idx),F(5,trl_idx,:)]           = utils.regain(f,E,toi,blperiod);
    
    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'symmetric','uniform','default','default');
    [R(6,trl_idx),F(6,trl_idx,:)]           = utils.regain(f,E,toi,blperiod);
    
    freqharm                                = [1:4]*tacsFreq;
    f                                       = artacs.dft.causal(r,freqharm,Fs,NumberPeriods);                  
    [R(7,trl_idx),F(7,trl_idx,:)]           = utils.regain(f,E,toi,blperiod);

    f                                       = artacs.template.stepwise(r,tacsFreq,Fs,'random');
    [R(8,trl_idx),F(8,trl_idx,:)]           = utils.regain(f,E,toi,blperiod);

    f                                       = clean_trials(datasample(1:size(clean_trials,1),1),:);
    [R(9,trl_idx),F(9,trl_idx,:)]           = utils.regain(f,E,toi,blperiod);

end


filter_names        = {'Causal Uniform','Causal Linear','Causal Exponential','Causal Gaussian','Symmetric Uniform','Symmetric Automatic','Causal Complex','Adaptive pPCA','Stim-free Bootstrap'};
Efun                = @(x,prm)mean(x,prm);

close all
for fidx = 1 : size(F,1)    
    figure
    set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
    hold on
    
    [h,p,ci] = ttest(squeeze(F(end,:,toi)));                
    h3 = patch([1:length(ci),length(ci):-1:1],[ci(1,:),fliplr(ci(2,:))],ones(1,length(ci)*2),'facecolor',red,'facealpha',0.25,'edgecolor',red,'edgealpha',0.1);                
    h1 = plot(mean(ci),'linewidth',1,'color',red);
    if fidx ~= size(F,1)        
        h2 = plot(squeeze(Efun(F(fidx,:,toi),2))','linewidth',2,'color',blue);
    end
    set(gca,'xlim',[1 length(toi)],'xtick',[1:50:4000],'xticklabel',[(-length(toi)/2+1):50:length(toi)/2])
    set(gca,'ylim',[-20 20],'ytick',[-1000:10:1000])
    title(filter_names{fidx})   
    print(gcf,[printfolder,'eva\ecg_td_',filter_names{fidx},'.png'],'-dpng')
end


%%%%%
clc
smpl_num    = 100;
koi = 0:0.01:1;
KxR = [];
KxA = [];
trl_num = size(F,2);
MR = [];
for fidx = 1 : size(R,1)
    
    kx              = ksdensity(R(fidx,:),koi,'width',0.05);
    kx              = kx./max(kx);   
    kx(kx<0.001)    = NaN;  
    KxR             = cat(1,KxR,kx);        
    
    bootstat    = [];    
    for bot_rep = 1 : 200
        bootstat = cat(2,bootstat,squeeze(Efun(F(fidx,datasample(1:trl_num,smpl_num),toi),2)));
    end
    
    r   = diag(corr(bootstat,repmat(E(1,toi),100,1)'));    
    MR = [MR,mean(r.^2)];
    kx  = ksdensity(r.^2,koi,'width',0.05);
    kx  = kx./max(kx);   
    kx(kx<0.001) =NaN;  
    KxA = cat(1,KxA,kx);       
end

%
close all
figure
set(gcf,'Position',[500 100 900 300],'paperpositionmode','auto')
for fidx = 1 : size(R,1)
    line(fidx+(.25*-KxR(fidx,:))-.2,koi,'color','k');
    line(fidx+(.25*+KxR(fidx,:))-.2,koi,'color','k');        
    m   = mean(R(fidx,:)); 
    hs  = line([fidx-.3,fidx-.1],[m,m],'color',blue,'linewidth',2);
    
    line(fidx+(.25*-KxA(fidx,:))+.2,koi,'color','k');
    line(fidx+(.25*+KxA(fidx,:))+.2,koi,'color','k');
    m   = MR(fidx);
    ha  = line([fidx+.1,fidx+.3],[m,m],'color',red,'linewidth',2);
    
end
set(gca,'ylim',[0 1],'ytick',[0:0.2:1],'yticklabel',[0:.2:1])
set(gca,'xlim',[0.5 fidx+.5],'xtick',1:fidx,'xticklabel',filter_names,'xaxislocation','bottom')
set(gca,'XTickLabelRotation',45)
lh = legend([hs,ha],'Single Trial',['Averaged n= ',num2str(smpl_num)]);
set(lh,'position',[0.1 .85 0 .130])
ylabel('Recovery (R²)')
grid on
print(gcf,[printfolder,'eva\ecg_R2.png'],'-dpng')
%% SIM PLOTS
setup           = generate.generic();  
[r,e,t]         = generate.recording(setup);

close all
figure
set(gcf,'Position',[100 100 500 400],'paperpositionmode','auto')
subplot(2,1,1)
hold on
h1 = plot(r,'color',[.8 .8 .8],'linewidth',2);
h2 = plot(e(1,:),'color',red);
set(gca,'xlim',[3751,4251],'xtick',[0:50:8000],'xticklabel',[-4000:50:4000])
set(gca,'ylim',[-21000 21000],'ytick',[-20000:10000:20000],'yticklabel',[-20:10:20])
ylabel('mV')
subplot(2,1,2)
hold on
h1 = plot(r,'color',[.8 .8 .8],'linewidth',2);
h2 = plot(e(1,:),'color',red);
set(gca,'xlim',[3751,4251],'xtick',[0:50:8000],'xticklabel',[-4000:50:4000])
set(gca,'ylim',[-3 3],'ytick',[-4:4])
ylabel('µV')
xlabel('ms')
lh = legend([h1,h2],'Artifacted Simulation','Clean Simulation');
set(lh,'position',[0.8 .85 .08 .08])
print(gcf,['.\repo\img\exemplary_generic.png'],'-dpng')


figure
set(gcf,'Position',[100 100 500 200],'paperpositionmode','auto')
hold on
h1 = plot(r,'color',[.8 .8 .8],'linewidth',2);
h2 = plot(e(1,:),'color',red);
set(gca,'xlim',[3750,4251],'xtick',[0:50:8000],'xticklabel',[-4000:50:4000])
set(gca,'ylim',[-21000 21000],'ytick',[-20000:10000:20000],'yticklabel',[-20:10:20])
ylabel('mV')
xlabel('ms')
print(gcf,[printfolder,'eva\sim_raw_1'],'-dpng')

figure
set(gcf,'Position',[100 100 500 200],'paperpositionmode','auto')
hold on
h1 = plot(r,'color',[.8 .8 .8],'linewidth',2);
h2 = plot(e(1,:),'color',red);
set(gca,'xlim',[3750,4251],'xtick',[0:50:8000],'xticklabel',[-4000:50:4000])
set(gca,'ylim',[-3 3],'ytick',[-4:4])
ylabel('µV')
xlabel('ms')
lh = legend([h1,h2],'Artifacted Simulation','Clean Simulation');
set(lh,'position',[0.7 .7 0.15 .25],'fontsize',12)
print(gcf,[printfolder,'eva\sim_raw_2'],'-dpng')
%

R               = [];    
F               = [];
setup           = generate.generic();  
tacsFreq        = 10;
NumberPeriods   = 10;
[r,e]           = generate.recording(setup);
E               = e(1,:);

for trl_idx = 1 : 100
    r                                       = generate.recording('generic');
    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'causal','uniform','default','default');
    [R(1,trl_idx),F(1,trl_idx,:)]           = utils.regain(f,E,toi,blperiod);

    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'causal','linear','default','default');
    [R(2,trl_idx),F(2,trl_idx,:)]           = utils.regain(f,E,toi,blperiod);

    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'causal','exp','default','default');
    [R(3,trl_idx),F(3,trl_idx,:)]           = utils.regain(f,E,toi,blperiod);

    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'causal','gauss','default','default');
    [R(4,trl_idx),F(4,trl_idx,:)]           = utils.regain(f,E,toi,blperiod);

    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'symmetric','uniform','default','default');
    [R(5,trl_idx),F(5,trl_idx,:)]           = utils.regain(f,E,toi,blperiod);

    freqharm                                = [1:4]*tacsFreq;
    f                                       = artacs.dft.causal(r,freqharm,Fs,NumberPeriods);                  
    [R(6,trl_idx),F(6,trl_idx,:)]           = utils.regain(f,E,toi,blperiod);

    f                                       = artacs.template.stepwise(r,tacsFreq,Fs,'random');
    [R(7,trl_idx),F(7,trl_idx,:)]           = utils.regain(f,E,toi,blperiod);
end



filter_names        = {'Causal Uniform','Causal Linear','Causal Exponential','Causal Gaussian','Symmetric Uniform','Causal DFT','Adaptive PCA','Stim-free Bootstrap'};
Efun                = @(x,prm)mean(x,prm);

close all
for fidx = 1 : size(F,1)    
    figure
    set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
    hold on
    
    h1 = plot(E(toi),'linewidth',1,'color',red);
    h2 = plot(squeeze(Efun(F(fidx,:,toi),2))','linewidth',2,'color',blue);
    set(gca,'xlim',[1 length(toi)],'xtick',[1:50:4000],'xticklabel',[(-length(toi)/2+1):50:length(toi)/2])
    set(gca,'ylim',[-3 3],'ytick',[-3:1:3])
    title(filter_names{fidx})   
    print(gcf,[printfolder,'eva\sim_td_',filter_names{fidx},'.png'],'-dpng')
end


%%%%%
clc
smpl_num    = 100;
koi = 0:0.01:1;
KxR = [];
KxA = [];
trl_num = size(F,2);
MR = [];
for fidx = 1 : size(R,1)
    
    kx              = ksdensity(R(fidx,:),koi,'width',0.05);
    kx              = kx./max(kx);   
    kx(kx<0.001)    = NaN;  
    KxR             = cat(1,KxR,kx);        
    
    bootstat    = [];    
    for bot_rep = 1 : 200
        bootstat = cat(2,bootstat,squeeze(Efun(F(fidx,datasample(1:trl_num,smpl_num),toi),2)));
    end
    
    r   = diag(corr(bootstat,repmat(E(1,toi),100,1)'));    
    MR = [MR,mean(r.^2)];
    kx  = ksdensity(r.^2,koi,'width',0.05);
    kx  = kx./max(kx);   
    kx(kx<0.001) =NaN;  
    KxA = cat(1,KxA,kx);       
end

%

figure
set(gcf,'Position',[500 100 900 300],'paperpositionmode','auto')
for fidx = 1 : size(R,1)
    line(fidx+(.25*-KxR(fidx,:))-.2,koi,'color','k');
    line(fidx+(.25*+KxR(fidx,:))-.2,koi,'color','k');        
    m   = mean(R(fidx,:)); 
    hs  = line([fidx-.3,fidx-.1],[m,m],'color',blue,'linewidth',2);
    
    line(fidx+(.25*-KxA(fidx,:))+.2,koi,'color','k');
    line(fidx+(.25*+KxA(fidx,:))+.2,koi,'color','k');
    m   = MR(fidx);
    ha  = line([fidx+.1,fidx+.3],[m,m],'color',red,'linewidth',2);
    
end
set(gca,'ylim',[0 1],'ytick',[0:0.2:1],'yticklabel',[0:.2:1])
set(gca,'xlim',[0.5 fidx+.5],'xtick',1:fidx,'xticklabel',filter_names,'xaxislocation','bottom')
set(gca,'XTickLabelRotation',45)
lh = legend([hs,ha],'Single Trial',['Averaged n= ',num2str(smpl_num)]);
set(lh,'position',[0.1 .85 0 .130])
ylabel('Recovery (R²)')
grid on
print(gcf,[printfolder,'eva\sim_R2.png'],'-dpng')

%% Frequency Response of Examples
close all

NumberPeriods   = 10;
tacsFreq        = 10;
Fs              = 1000;

artacs.kernel.response(artacs.kernel.causal(NumberPeriods,tacsFreq,Fs,'ave'),1)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
set(gca,'ylim',[-.75 1.1])
print(gcf,[printfolder,'causal/kernel_ave.png'],'-dpng')

artacs.kernel.response(artacs.kernel.causal(NumberPeriods,tacsFreq,Fs,'linear'),1)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
set(gca,'ylim',[-.75 1.1])
print(gcf,[printfolder,'causal/kernel_linear.png'],'-dpng')

artacs.kernel.response(artacs.kernel.causal(NumberPeriods,tacsFreq,Fs,'exp'),1)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
set(gca,'ylim',[-.75 1.1])
print(gcf,[printfolder,'causal/kernel_exp.png'],'-dpng')

artacs.kernel.response(artacs.kernel.causal(NumberPeriods,tacsFreq,Fs,'gauss'),1)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
set(gca,'ylim',[-.75 1.1])
print(gcf,[printfolder,'causal/kernel_gauss.png'],'-dpng')
% --------------------

artacs.kernel.response(artacs.kernel.causal(NumberPeriods,tacsFreq,Fs,'ave'),2,30)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
print(gcf,[printfolder,'causal/mag_ave.png'],'-dpng')

artacs.kernel.response(artacs.kernel.causal(NumberPeriods,tacsFreq,Fs,'linear'),2,30)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
print(gcf,[printfolder,'causal/mag_linear.png'],'-dpng')

artacs.kernel.response(artacs.kernel.causal(NumberPeriods,tacsFreq,Fs,'exp'),2,30)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
print(gcf,[printfolder,'causal/mag_exp.png'],'-dpng')

artacs.kernel.response(artacs.kernel.causal(NumberPeriods,tacsFreq,Fs,'gauss'),2,30)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
print(gcf,[printfolder,'causal/mag_gauss.png'],'-dpng')

%%

artacs.kernel.response(artacs.kernel.symmetric(NumberPeriods,tacsFreq,Fs,'ave'),1)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
set(gca,'ylim',[-.75 1.1])
print(gcf,[printfolder,'sym/kernel_ave.png'],'-dpng')

artacs.kernel.response(artacs.kernel.symmetric(NumberPeriods,tacsFreq,Fs,'linear'),1)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
set(gca,'ylim',[-.75 1.1])
print(gcf,[printfolder,'sym/kernel_linear.png'],'-dpng')

artacs.kernel.response(artacs.kernel.symmetric(NumberPeriods,tacsFreq,Fs,'exp'),1)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
set(gca,'ylim',[-.75 1.1])
print(gcf,[printfolder,'sym/kernel_exp.png'],'-dpng')

artacs.kernel.response(artacs.kernel.symmetric(NumberPeriods,tacsFreq,Fs,'gauss'),1)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
set(gca,'ylim',[-.75 1.1])
print(gcf,[printfolder,'sym/kernel_gauss.png'],'-dpng')

% --------------------

artacs.kernel.response(artacs.kernel.symmetric(NumberPeriods,tacsFreq,Fs,'ave'),2,30)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
print(gcf,[printfolder,'sym/mag_ave.png'],'-dpng')

artacs.kernel.response(artacs.kernel.symmetric(NumberPeriods,tacsFreq,Fs,'linear'),2,30)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
print(gcf,[printfolder,'sym/mag_linear.png'],'-dpng')

artacs.kernel.response(artacs.kernel.symmetric(NumberPeriods,tacsFreq,Fs,'exp'),2,30)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
print(gcf,[printfolder,'sym/mag_exp.png'],'-dpng')

artacs.kernel.response(artacs.kernel.symmetric(NumberPeriods,tacsFreq,Fs,'gauss'),2,30)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
print(gcf,[printfolder,'sym/mag_gauss.png'],'-dpng')

%%

artacs.kernel.response(artacs.kernel.create(NumberPeriods,tacsFreq,Fs,'symmetric','uniform','default','inc',1),1)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
set(gca,'ylim',[-.75 1.1])
print(gcf,[printfolder,'inv/kernel_ave.png'],'-dpng')

artacs.kernel.response(artacs.kernel.create(NumberPeriods,tacsFreq,Fs,'symmetric','linear','default','inc',1),1)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
set(gca,'ylim',[-.75 1.1])
print(gcf,[printfolder,'inv/kernel_linear.png'],'-dpng')

artacs.kernel.response(artacs.kernel.create(NumberPeriods,tacsFreq,Fs,'symmetric','exp','default','inc',1),1)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
set(gca,'ylim',[-.75 1.1])
print(gcf,[printfolder,'inv/kernel_exp.png'],'-dpng')

artacs.kernel.response(artacs.kernel.create(NumberPeriods,tacsFreq,Fs,'symmetric','gauss','default','inc',1),1)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
set(gca,'ylim',[-.75 1.1])
print(gcf,[printfolder,'inv/kernel_gauss.png'],'-dpng')

% --------------------

artacs.kernel.response(artacs.kernel.create(NumberPeriods,tacsFreq,Fs,'symmetric','uniform','default','inc',1),2,30)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
print(gcf,[printfolder,'inv/mag_ave.png'],'-dpng')

artacs.kernel.response(artacs.kernel.create(NumberPeriods,tacsFreq,Fs,'symmetric','linear','default','inc',1),2,30)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
print(gcf,[printfolder,'inv/mag_linear.png'],'-dpng')

artacs.kernel.response(artacs.kernel.create(NumberPeriods,tacsFreq,Fs,'symmetric','exp','default','inc',1),2,30)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
print(gcf,[printfolder,'inv/mag_exp.png'],'-dpng')

artacs.kernel.response(artacs.kernel.create(NumberPeriods,tacsFreq,Fs,'symmetric','gauss','default','inc',1),2,30)
set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
print(gcf,[printfolder,'inv/mag_gauss.png'],'-dpng')

%% Kernel Tuning
tau     = [0,5,100,1000];
close all

for tau_idx = 1 : length(tau)
    figure
    hold on
    set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
    artacs.kernel.response(artacs.kernel.causal(NumberPeriods,tacsFreq,Fs,'exp',tau(tau_idx)),1,[],gcf)    
    set(gca,'ylim',[-1.1 1.1])
    %set(get(gca,'children'),'Markerfacecolor',cmap{tau_idx});    
    print(gcf,[printfolder,'tau/kernel_exp_',num2str(tau(tau_idx)),'.png'],'-dpng')
    
    figure
    hold on
    set(gcf,'Position',[100 500 400 300],'paperpositionmode','auto')
    artacs.kernel.response(artacs.kernel.causal(NumberPeriods,tacsFreq,Fs,'exp',tau(tau_idx)),2,30,gcf)    
   % set(get(gca,'children'),'Color',cmap{tau_idx});    
    print(gcf,[printfolder,'tau/mag_exp_',num2str(tau(tau_idx)),'.png'],'-dpng')
end

%%
close all
NumberPeriods = [1,10,50];
for N = 1 : length(NumberPeriods)
    figure
    hold on
    set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
    artacs.kernel.response(artacs.kernel.causal(NumberPeriods(N),tacsFreq,Fs,'exp'),1,[],gcf)    
    set(gca,'ylim',[-1.1 1.1])
    %set(get(gca,'children'),'Markerfacecolor',cmap{tau_idx});    
    print(gcf,[printfolder,'np/kernel_exp_',num2str(NumberPeriods(N)),'.png'],'-dpng')
    
    figure
    hold on
    set(gcf,'Position',[100 500 400 300],'paperpositionmode','auto')
    artacs.kernel.response(artacs.kernel.causal(NumberPeriods(N),tacsFreq,Fs,'exp'),2,30,gcf)    
   % set(get(gca,'children'),'Color',cmap{tau_idx});    
    print(gcf,[printfolder,'np/mag_exp_',num2str(NumberPeriods(N)),'.png'],'-dpng')
           
end


