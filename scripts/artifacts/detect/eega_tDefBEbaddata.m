% -------------------------------------------------------------------------
% This functions defines which epoch are bad based on the amount of
% rejected data
%
% INPUTS
% EEG   EEG structure
% BadData   structure speficfing the maximun amount of bad data
%   - limBCTa   maximun data rejected in BCT
%   - limBTa    maximun data rejected in BT
%   - limBCa    maximun data rejected in BC
%   - limCCTa   maximun data corrected in CCT
%   - limBCTr   maximun data rejected in BCT (relative limit)
%   - limBTr    maximun data rejected in BT (relative limit)
%   - limBCr    maximun data rejected in BC (relative limit)
%   - limCCTr   maximun data corrected in CCT (relative limit)
%
% OPTIONAL INPUTS
%   - keeppre       keep previuos values
%   - where         time limits where to look for bad data
%   - maxloops      maximun numbers of loops
%   - plot          plot the rejection
%   - savefigure    save the figure with the rejection
%
% OUTPUTS
%   EEG     output data
%   BE     bad epochs, logical indexes 
%
% -------------------------------------------------------------------------

function [ EEG, BE ] = eega_tDefBEbaddata( EEG, BadData, varargin )

fprintf('### Identifying Bad Epochs ###\n' )

%% ------------------------------------------------------------------------
%% Parameters
P.keeppre       = 1;
P.log           = 0;
P.maxloops      = 1;
P.plot          = 0;
P.savefigure    = 0;
P.where         = [];

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tDefBEart: Non recognized inputs')
end

if ~isfield(BadData,'limBCTr')
    BadData.limBCTr = [];
end
if ~isfield(BadData,'limBTr')
    BadData.limBTr = [];
end
if ~isfield(BadData,'limBCr')
    BadData.limBCr = [];
end
if ~isfield(BadData,'limCCTr')
    BadData.limCCTr = [];
end
limRel = {BadData.limBCTr, BadData.limBTr, BadData.limBCr, BadData.limCCTr};
limAbs = {BadData.limBCTa, BadData.limBTa, BadData.limBCa, BadData.limCCTa};


%% ------------------------------------------------------------------------
%% Identify bad epochs
nEl  = size(EEG.data,1);
nS   = size(EEG.data,2);
nEp  = size(EEG.data,3);

if P.keeppre && isfield(EEG.artifacts,'BE')
    BEold = EEG.artifacts.BE;
else
    BEold = false(1,1,nEp);
end
if ~isfield(EEG.artifacts,'BEm')
    EEG.artifacts.BEm = false(1,1,nEp);
end
if ~isfield(EEG.artifacts,'BCT')
    EEG.artifacts.BCT = false(nEl,nS,nEp);
end
if ~isfield(EEG.artifacts,'BT')
    EEG.artifacts.BT = false(1,nS,nEp);
end
if ~isfield(EEG.artifacts,'BC')
    EEG.artifacts.BC = false(nEl,1,nEp);
end
if ~isfield(EEG.artifacts,'CCT')
    EEG.artifacts.CCT = false(nEl,nS,nEp);
end

% find the times to consider
if isempty(P.where)
    P.where = [EEG.times(1) EEG.times(end)];
end
idxtime = EEG.times>=P.where(1) & EEG.times<=P.where(2);
nS = sum(idxtime);

% find bad epochs
% ElBadAll = all(EEG.artifacts.BC,3);
% nElBadAll = sum(ElBadAll);
% nElAct = nEl - nElBadAll;
% TOT = nan(nEp,4);
% TOT(:,1) = squeeze(sum(sum(EEG.artifacts.BCT(~ElBadAll,idxtime,:),1),2) / (nS*nElAct));
% TOT(:,2) = squeeze(sum(EEG.artifacts.BT(:,idxtime,:),2) / nS);
% TOT(:,3) = squeeze(sum(EEG.artifacts.BC(~ElBadAll,:,:),1) / nElAct);
% TOT(:,4) = squeeze(sum(sum(EEG.artifacts.CCT(~ElBadAll,idxtime,:),1),2) / (nS*nElAct));
TOT = nan(nEp,4);
TOT(:,1) = squeeze(sum(sum(EEG.artifacts.BCT(:,idxtime,:),1),2) / (nS*nEl));
TOT(:,2) = squeeze(sum(EEG.artifacts.BT(:,idxtime,:),2) / nS);
TOT(:,3) = squeeze(sum(EEG.artifacts.BC(:,:,:),1) / nEl);
TOT(:,4) = squeeze(sum(sum(EEG.artifacts.CCT(:,idxtime,:),1),2) / (nS*nEl));
if P.log
%     TOT(TOT(:,1)==0,1) = 1/ (nS*nElAct);
%     TOT(TOT(:,2)==0,2) = 1/ nS;
%     TOT(TOT(:,3)==0,3) = 1/ nElAct;
%     TOT(TOT(:,4)==0,4) = 1/ (nS*nElAct);
    TOT(TOT(:,1)==0,1) = 1/ (nS*nEl);
    TOT(TOT(:,2)==0,2) = 1/ nS;
    TOT(TOT(:,3)==0,3) = 1/ nEl;
    TOT(TOT(:,4)==0,4) = 1/ (nS*nEl);
    TOT = log(TOT);
    for i=1:4
        limAbs{i} = log(limAbs{i});
    end
end

BE = BEold(:) | EEG.artifacts.BEm(:);
BEd = false(size(BE));
ok=0;
ci=1;
while ~ok && ci<=P.maxloops
    thresh = ones(1,4);
    for i=1:4
        if ~isempty(limRel{i})
            P75 = prctile(TOT(~BE,i),75);
            P25 = prctile(TOT(~BE,i),25);
            thresh(i) = P75 + limRel{i} .* (P75-P25);
            if ~isempty(limAbs{i}) && length(limAbs{i})==2
                if thresh(i)<limAbs{i}(1); thresh(i) = limAbs{i}(1);
                elseif thresh(i)>limAbs{i}(2); thresh(i) = limAbs{i}(2);
                end
            end
        else
            thresh(i) = limAbs{i}(1);
        end
    end
    R = TOT > repmat(thresh,[nEp 1]);
    
    if all( (any(R,2) | BE ) == BE) % check if new data was rejected
        ok=1;
    end
    
    BEd = ( BEd | any(R,2) );
    BE  = ( BE  | any(R,2) );
    
    ci=ci+1;
end

%% ------------------------------------------------------------------------
%% Display rejected data
BEnew = (BE(:) & ~BEold(:) & ~EEG.artifacts.BEm(:));

fprintf('--> Rejected epochs by this algorithm: %03d out of %d : (%5.1f%% ) %s\n', sum(BEd), nEp, sum(BEd)/nEp*100, num2str(find(BEd(:)')) )
fprintf('       - BCT threshold %4.3f, trials %d (%5.2f%%): %s\n', thresh(1), sum(R(:,1)), sum(R(:,1))/nEp*100, num2str(find(R(:,1)')) )
fprintf('       - BT  threshold %4.3f, trials %d (%5.2f%%): %s\n', thresh(2), sum(R(:,2)), sum(R(:,2))/nEp*100, num2str(find(R(:,2)')) )
fprintf('       - BC  threshold %4.3f, trials %d (%5.2f%%): %s\n', thresh(3), sum(R(:,3)), sum(R(:,3))/nEp*100, num2str(find(R(:,3)')) )
fprintf('       - CCT threshold %4.3f, trials %d (%5.2f%%): %s\n', thresh(4), sum(R(:,4)), sum(R(:,4))/nEp*100, num2str(find(R(:,4)')) )

fprintf('--> Total rejected epochs:             %03d out of %d (%5.1f%% ) %s\n', sum(BE), nEp, sum(BE)/nEp*100, num2str(find(BE(:)')) )
fprintf('--> New rejected epochs:               %03d out of %d (%5.1f%% ) %s\n', sum(BEnew), nEp, sum(BEnew)/nEp*100, num2str(find(BEnew(:)')) )


fprintf('\n')

%% ------------------------------------------------------------------------
%% Update the rejection matrix
EEG.artifacts.BE = permute(BE,[3 2 1]);
EEG.reject.rejmanual = permute(EEG.artifacts.BE,[1 3 2]);
if exist('eega_summarypp','file')==2
    EEG = eega_summarypp(EEG);
end

%% ------------------------------------------------------------------------
%% Plot
if P.plot
    plottrialsrej(BE,BEd,BEold,EEG.artifacts,TOT,thresh,P,EEG.filename,EEG.filepath)
end

end

%% ------------------------------------------------------------------------

function plottrialsrej(BE,BEd,BEold,Art,TOT,thresh,P,filename,filepath)

D=Art.BCT;
D(repmat(Art.BC,[1 size(D,2) 1]))=1;
D(repmat(Art.BT,[size(D,1) 1 1]))=1;
D = permute(sum(D,1),[2 3 1]);

AXxLim = [0 0.05 0.67 0.7 0.98 1];
AXyLim = [0 0.1 0.90 1];
axbox  = 0.055;
axboxm  = 0.015;

ppp ={'BCT' 'BT' 'BC' 'CCT'};

col_good    = [0.1953    0.8008    0.1953];
col_new     = [1.0000    0.8398         0];
col_old     = [0.2539    0.4102    0.8789];
col_both    = [0.8594    0.0781    0.2344];

Egood   = ~BE(:);
Enew    = BEd & ~BEold(:);
Eold    = BEold(:) & ~BEd(:);
Eboth   = BEold(:) & BEd(:);

EEEE = [Egood Enew Eold Eboth];

figure('Position',[100 100 1200 500])

% trial distance
% ---------------
axes('Position', [AXxLim(2) AXyLim(2) AXxLim(3)-AXxLim(2) AXyLim(3)-AXyLim(2)])

% yLim = [0 ( prctile(D(:),75) + 2.5 * (prctile(D(:),75) - prctile(D(:),25)) )];
yLim = [0 size(D,2)*0.35];

triaslorder=[];
L={};
C=[];
Nlim=0;
i=0;
if any(EEEE(:,1))
    i=i+1;
    triaslorder = cat(1,triaslorder,find(EEEE(:,1)));
    Nlim(i+1) = Nlim(i)+sum(EEEE(:,1));
    L{i} = 'good epoch';
    C(i,:)= col_good;
end
if any(EEEE(:,2))
    i=i+1;
    triaslorder = cat(1,triaslorder,find(EEEE(:,2)));
    Nlim(i+1) = Nlim(i)+sum(EEEE(:,2));
    L{i} = 'rejected new';
    C(i,:)= col_new;
end
if any(EEEE(:,3))
    i=i+1;
    triaslorder = cat(1,triaslorder,find(EEEE(:,3)));
    Nlim(i+1) = Nlim(i)+sum(EEEE(:,3));
    L{i} = 'rejected old';
    C(i,:)= col_old;
end
if any(EEEE(:,4))
    i=i+1;
    triaslorder = cat(1,triaslorder,find(EEEE(:,4)));
    Nlim(i+1) = Nlim(i)+sum(EEEE(:,4));
    L{i} = 'rejected both';
    C(i,:)= col_both;
end
Nlim=Nlim+0.5;

Dorder = D(:,triaslorder);
imagesc(Dorder')
colormap(jet)
caxis(yLim)
set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[])
xlabel('time')
ylabel('trial')
n=-round(size(Dorder,1)*0.02);
xlim([n size(Dorder,1)])
H=[];
for j=2:length(Nlim)
    v=[n Nlim(j); 0.5 Nlim(j); 0.5 Nlim(j-1); n Nlim(j-1)];
    h=patch(v(:,1),v(:,2),C(j-1,:),'EdgeColor','none');
    H{j-1}=h(1);
end
legend([H{:}]',L{:},'Location','southoutside','Orientation','horizontal')
title('bad data per epoch')
colorbar

% boxplots
% --------
X0=AXxLim(4);
for i=1:size(TOT,2)
    axes('Position', [X0 AXyLim(2) axbox AXyLim(3)-AXyLim(2)])
    hold on
    di=TOT(:,i);
    digood = di;
    digood(BE)=nan;
    boxplot([di digood],'Labels',{'all' 'good'},'PlotStyle','compact', 'Colors' ,[0 0 0])
    set(gca,'XTickLabelRotation',90)
%     set(gca,'YTickLabel',[])
    line(xlim, [thresh(i) thresh(i)],'color',[1 0 0])
    title(ppp{i})
    X0=X0+(axbox+axboxm);
end

if P.plot && P.savefigure
    [~,figurename,~]=fileparts(filename);
    figurename = ['rejBEart_' figurename ];
    savefig(fullfile(filepath,figurename))
end
end

