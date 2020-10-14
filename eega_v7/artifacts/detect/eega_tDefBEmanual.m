function [ EEG, BE ] = eega_tDefBEmanual( EEG, varargin )

fprintf('### Manual rejection of epochs ###\n' )

%% ------------------------------------------------------------------------
%% Parameters
P.keeppre       = 0;
P.plot          = 1;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tDefBTBC: Non recognized inputs')
end

%% ------------------------------------------------------------------------
%% Identify bad epochs
nEl  = size(EEG.data,1);
nS   = size(EEG.data,2);
nEp  = size(EEG.data,3);

if isfield(EEG.artifacts,'BEm')
    BEmold = EEG.artifacts.BEm;
else
    BEmold = false(nEp,1);
end
if isfield(EEG.artifacts,'BE')
    BEold = EEG.artifacts.BE;
else
    BEold = false(nEp,1);
end

if isfield(EEG.artifacts,'BCT')
    BCT = EEG.artifacts.BCT;
else
	BCT = false(nEl,nS,nEp);
end
if isfield(EEG.artifacts,'BT')
	BT = EEG.artifacts.BT;
else
    BT = false(1,nS,nEp);
end
if isfield(EEG.artifacts,'BC')
	BC = EEG.artifacts.BC;
else
	BC = false(nEl,1,nEp);
end
if isfield(EEG.artifacts,'CCT')
	CCT = EEG.artifacts.CCT;
else
	CCT = false(nEl,nS,nEp);
end

BEm = false(nEp,1);
BE = BEold;
if P.plot
    BEm = false(nEp,1);
    ax_bct = [600 300];
    ax_bc = [30 300];
    ax_bt = [600 30];
    ax_be = [600 50];
    ax_erp = [600 300];
    fig_ = [750 1300];
    ax_bct = ax_bct./fig_;
    ax_bc = ax_bc./fig_;
    ax_bt = ax_bt./fig_;
    ax_erp = ax_erp./fig_;
    ax_be = ax_be./fig_;
    nnn=11;
    
    cmap_bct = [0 0 0; 0.75 0.25 0];
    cmap_cct = [0 0 0; 0 0.75 0.25];
    
    d=EEG.data;
    d(BCT) = nan;
    yl = 4*(prctile(d(:),75)-prctile(d(:),25));
    yl = [-yl yl];
        
    for ep=1:nEp
        
        hfep = figure('Position',[1 1 fig_]);
        
        hax_be = axes('Position',[[50 1200]./fig_ ax_be]);
        text(0,0.8,sprintf('Epoch %d',ep),'FontSize',16)
        text(0,0.2,sprintf('BE = %d',BEold(ep)),'FontSize',16)
        set(hax_be,'Visible','off')
        
        hax_bct = axes('Position',[[50 850]./fig_ ax_bct]);
        imagesc(BCT(:,:,ep)), colormap(hax_bct,cmap_bct)
        set(gca,'XTick',[])
        ylabel('channel')
        
        hax_bt = axes('Position',[[50 800]./fig_ ax_bt]);
        imagesc(BT(:,:,ep)), colormap(hax_bt,cmap_bct)
        set(gca,'XTick',1:floor(nS/nnn):nS)
        set(gca,'XTickLabel',EEG.times(1:floor(nS/nnn):nS))
        xlabel('time')
        
        hax_bc = axes('Position',[[670 850]./fig_ ax_bc]);
        imagesc(BC(:,:,ep)), colormap(hax_bc,cmap_bct)
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        
        hax_cct = axes('Position',[[50 450]./fig_ ax_bct]);
        imagesc(CCT(:,:,ep)), colormap(hax_cct,cmap_cct)
        ylabel('channel')
        set(gca,'XTick',1:floor(nS/nnn):nS)
        set(gca,'XTickLabel',EEG.times(1:floor(nS/nnn):nS))
        
        hax_erp = axes('Position',[[50 50]./fig_ ax_erp]);
        hold on
        plot(EEG.times,EEG.data(:,:,ep)','Color',[0 0 0])
        d=EEG.data(:,:,ep);
        d(~CCT(:,:,ep)) = nan;
        plot(EEG.times,d,'Color',[0 0.75 0.25])
        d=EEG.data(:,:,ep);
        d(~BCT(:,:,ep)) = nan;
        plot(EEG.times,d','Color',[0.75 0.25 0])
        d=EEG.data(:,:,ep);
        d(~(BCT(:,:,ep) & ~CCT(:,:,ep)) ) = nan;
        plot(EEG.times,d','Color',[0.05 0.05 0.95])
        ylim(yl)
        
        ok = 0;
        while ~ok
            rej = input(sprintf('Do you want to reject epoch %d? 0=nothing, 1=reject, 2=include: ',ep));
            if rej==0
                ok=1;
            elseif rej==1
                ok=1;
                BEm(ep) = 1;
            elseif rej==2
                ok=1;
                BE(ep) = 0;
            end
        end
        close(hfep)
    end
else
    badep = input('Epochs to reject manually:');
    BEm(badep)=true;
end

EEG.artifacts.BEm = BEm;
EEG.artifacts.BE = BEm | BE;

