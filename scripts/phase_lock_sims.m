%% OSCI SIMULATION _phaselocked

setup           = generate.generic();  
setup.erpMagnitude = 0;

setup.tacsPhase = 0;
%setup.tacsPhase = 'random';
setup.eoPhase   = 0;
tacsFreq        = 10;
NumberPeriods   = 10;
[r,e,t]         = generate.recording(setup);

close all
figure
set(gcf,'Position',[100 100 500 200],'paperpositionmode','auto')
hold on
h1 = plot(r,'color',[.8 .8 .8],'linewidth',2);
h2 = plot(e(2,:),'color',red);
set(gca,'xlim',[3750,4251],'xtick',[0:50:8000],'xticklabel',[-4000:50:4000])
set(gca,'ylim',[-21000 21000],'ytick',[-20000:10000:20000],'yticklabel',[-20:10:20])
ylabel('mV')
xlabel('ms')
print(gcf,[printfolder,'eva\eo_raw_1'],'-dpng')

figure
set(gcf,'Position',[100 100 500 200],'paperpositionmode','auto')
hold on
h1 = plot(r,'color',[.8 .8 .8],'linewidth',2);
h2 = plot(e(2,:),'color',red);
set(gca,'xlim',[3750,4251],'xtick',[0:50:8000],'xticklabel',[-4000:50:4000])
set(gca,'ylim',[-3 3],'ytick',[-4:4])
ylabel('µV')
xlabel('ms')
lh = legend([h1,h2],'Artifacted Simulation','Clean Simulation');
set(lh,'position',[0.7 .7 0.15 .25],'fontsize',12)
print(gcf,[printfolder,'eva\eo_raw_2'],'-dpng')
%

R               = [];    
F               = [];
[r,e]           = generate.recording(setup);
E               = abs(hilbert(e(2,:)));

for trl_idx = 1 : 100
    r                                       = generate.recording(setup);
    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'causal','uniform','default','default');
    [R(1,trl_idx),F(1,trl_idx,:)]           = utils.regain(abs(hilbert(f)),E,toi,blperiod);

    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'causal','linear','default','default');
    [R(2,trl_idx),F(2,trl_idx,:)]           = utils.regain(abs(hilbert(f)),E,toi,blperiod);

    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'causal','exp','default','default');
    [R(3,trl_idx),F(3,trl_idx,:)]           = utils.regain(abs(hilbert(f)),E,toi,blperiod);

    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'causal','gauss','default','default');
    [R(4,trl_idx),F(4,trl_idx,:)]           = utils.regain(abs(hilbert(f)),E,toi,blperiod);

    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'symmetric','uniform','default','default');
    [R(5,trl_idx),F(5,trl_idx,:)]           = utils.regain(abs(hilbert(f)),E,toi,blperiod);

    freqharm                                = [1:4]*tacsFreq;
    f                                       = artacs.dft.causal(r,freqharm,Fs,NumberPeriods);                  
    [R(6,trl_idx),F(6,trl_idx,:)]           = utils.regain(abs(hilbert(f)),E,toi,blperiod);

    f                                       = artacs.template.stepwise(r,tacsFreq,Fs,'random');
    [R(7,trl_idx),F(7,trl_idx,:)]           = utils.regain(abs(hilbert(f)),E,toi,blperiod);
end



filter_names        = {'Causal Uniform','Causal Linear','Causal Exponential','Causal Gaussian','Symmetric Uniform','Causal DFT','Adaptive PCA','Stim-free Bootstrap'};
Efun                = @(x,prm)mean(x,prm);

close all
for fidx = 1 : size(F,1)    
    figure
    set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
    hold on
    
    h1 = plot(E(toi),'linewidth',1,'color',red);
    x = squeeze(Efun(abs(hilbert(F(fidx,:,toi))),2))';
    h2 = plot(x-mean(x),'linewidth',2,'color',blue);
    set(gca,'xlim',[1 length(toi)],'xtick',[1:50:4000],'xticklabel',[(-length(toi)/2+1):50:length(toi)/2])
    set(gca,'ylim',[-3 3],'ytick',[-3:1:3])
    title(filter_names{fidx})   
    print(gcf,[printfolder,'eva\eo_td_',filter_names{fidx},'.png'],'-dpng')
end


%%%%%
clc
smpl_num    = 100;
koi         = 0:0.01:1;
KxR         = [];
KxA         = [];
trl_num     = size(F,2);
MR          = [];
for fidx = 1 : size(R,1)    
    kx              = ksdensity(R(fidx,:),koi,'width',0.05);
    kx              = kx./max(kx);   
    kx(kx<0.001)    = NaN;  
    KxR             = cat(1,KxR,kx);        
    
    bootstat    = [];    
    for bot_rep = 1 : 200
        bootstat = cat(2,bootstat,squeeze(Efun(abs(hilbert(F(fidx,datasample(1:trl_num,smpl_num),toi))),2)));
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
print(gcf,[printfolder,'eva\eo_R2.png'],'-dpng')
%% OSCI SIMULATION _random phase

setup           = generate.generic();  
setup.erpMagnitude = 0;

setup.tacsPhase = 0;
setup.tacsPhase = 'random';
tacsFreq        = 10;
NumberPeriods   = 10;
[r,e,t]         = generate.recording(setup);

close all
figure
set(gcf,'Position',[100 100 500 200],'paperpositionmode','auto')
hold on
h1 = plot(r,'color',[.8 .8 .8],'linewidth',2);
h2 = plot(e(2,:),'color',red);
set(gca,'xlim',[3750,4251],'xtick',[0:50:8000],'xticklabel',[-4000:50:4000])
set(gca,'ylim',[-21000 21000],'ytick',[-20000:10000:20000],'yticklabel',[-20:10:20])
ylabel('mV')
xlabel('ms')
print(gcf,[printfolder,'eva\random_eo_raw_1'],'-dpng')

figure
set(gcf,'Position',[100 100 500 200],'paperpositionmode','auto')
hold on
h1 = plot(r,'color',[.8 .8 .8],'linewidth',2);
h2 = plot(e(2,:),'color',red);
set(gca,'xlim',[3750,4251],'xtick',[0:50:8000],'xticklabel',[-4000:50:4000])
set(gca,'ylim',[-3 3],'ytick',[-4:4])
ylabel('µV')
xlabel('ms')
lh = legend([h1,h2],'Artifacted Simulation','Clean Simulation');
set(lh,'position',[0.7 .7 0.15 .25],'fontsize',12)
print(gcf,[printfolder,'eva\random_eo_raw_2'],'-dpng')
%

R               = [];    
F               = [];
[r,e]           = generate.recording(setup);
E               = abs(hilbert(e(2,:)));

for trl_idx = 1 : 100
    r                                       = generate.recording(setup);
    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'causal','uniform','default','default');
    [R(1,trl_idx),F(1,trl_idx,:)]           = utils.regain(abs(hilbert(f)),E,toi,blperiod);

    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'causal','linear','default','default');
    [R(2,trl_idx),F(2,trl_idx,:)]           = utils.regain(abs(hilbert(f)),E,toi,blperiod);

    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'causal','exp','default','default');
    [R(3,trl_idx),F(3,trl_idx,:)]           = utils.regain(abs(hilbert(f)),E,toi,blperiod);

    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'causal','gauss','default','default');
    [R(4,trl_idx),F(4,trl_idx,:)]           = utils.regain(abs(hilbert(f)),E,toi,blperiod);

    f                                       = artacs.kernel.run(r,NumberPeriods,tacsFreq,Fs,'symmetric','uniform','default','default');
    [R(5,trl_idx),F(5,trl_idx,:)]           = utils.regain(abs(hilbert(f)),E,toi,blperiod);

    freqharm                                = [1:4]*tacsFreq;
    f                                       = artacs.dft.causal(r,freqharm,Fs,NumberPeriods);                  
    [R(6,trl_idx),F(6,trl_idx,:)]           = utils.regain(abs(hilbert(f)),E,toi,blperiod);

    f                                       = artacs.template.stepwise(r,tacsFreq,Fs,'random');
    [R(7,trl_idx),F(7,trl_idx,:)]           = utils.regain(abs(hilbert(f)),E,toi,blperiod);
end



filter_names        = {'Causal Uniform','Causal Linear','Causal Exponential','Causal Gaussian','Symmetric Uniform','Causal DFT','Adaptive PCA','Stim-free Bootstrap'};
Efun                = @(x,prm)mean(x,prm);

close all
for fidx = 1 : size(F,1)    
    figure
    set(gcf,'Position',[100 100 400 300],'paperpositionmode','auto')
    hold on
    
    h1 = plot(E(toi),'linewidth',1,'color',red);
    x = squeeze(Efun(abs(hilbert(F(fidx,:,toi))),2))';
    h2 = plot(x-mean(x),'linewidth',2,'color',blue);
    set(gca,'xlim',[1 length(toi)],'xtick',[1:50:4000],'xticklabel',[(-length(toi)/2+1):50:length(toi)/2])
    set(gca,'ylim',[-3 3],'ytick',[-3:1:3])
    title(filter_names{fidx})   
    print(gcf,[printfolder,'eva\random_eo_td_',filter_names{fidx},'.png'],'-dpng')
end


%%%%%
clc
smpl_num    = 100;
koi         = 0:0.01:1;
KxR         = [];
KxA         = [];
trl_num     = size(F,2);
MR          = [];
for fidx = 1 : size(R,1)    
    kx              = ksdensity(R(fidx,:),koi,'width',0.05);
    kx              = kx./max(kx);   
    kx(kx<0.001)    = NaN;  
    KxR             = cat(1,KxR,kx);        
    
    bootstat    = [];    
    for bot_rep = 1 : 200
        bootstat = cat(2,bootstat,squeeze(Efun(abs(hilbert(F(fidx,datasample(1:trl_num,smpl_num),toi))),2)));
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
print(gcf,[printfolder,'eva\random_eo_R2.png'],'-dpng')
%%