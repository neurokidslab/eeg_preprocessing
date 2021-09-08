% this function uses eegplotto plot the rejection matrixes
% By default it plots BT and BC(all recording)
% It can also plot BTC (bad times x channles)
% this can be defined by setting the optional inputs 
% - plotbct (defaoult 0)
% - plotbt (defaoult 1)
% - plotbc (defaoult 1)

function eega_plot_rejection(EEG, plotbct, plotbt, plotbc, winlength)

if nargin<5
    winlength = min(120, size(EEG.data,2)/EEG.srate);
end
if nargin<4
    plotbc = 1;
end
if nargin<3
    plotbt = 1;
end
if nargin<2
    plotbct = 0;
end

[nEl nSm nEp] = size(EEG.data);
if nEp>1
    error('This function is not prepared for epoched data')
end

% create the matrix holding what to mark as rejected
winrej = [];

if plotbct  % plot in red the rejected data (white background)
    bct = reshape(EEG.artifacts.BCT,[nEl nSm*nEp]);
    for el=1:nEl
        btel = bct(el,:);
        btel_ini = find([btel(1) diff(btel,1,2)==1])';
        btel_fin = find([diff(btel,1,2)==-1 btel(end)])';
        for i=1:length(btel_ini)
            ch = false(1, nEl);
            ch(el) = 1;
            winreji = [btel_ini(i) btel_fin(i) 1 1 1 ch];
            winrej = cat(1,winrej,winreji);
        end
    end
end

if plotbt % plot with orange background the rejected times 
    bt = reshape(EEG.artifacts.BT,[1 nSm*nEp]);
    bt_ini = find([bt(1) diff(bt,1,2)==1])';
    bt_fin = find([diff(bt,1,2)==-1 bt(end)])';
    for i=1:length(bt_ini)
        winreji = [bt_ini(i) bt_fin(i) 0.8, 0.50, 0.2 true(1, nEl)];
        winrej = cat(1,winrej,winreji);
    end
end

if plotbc % plot in red the rejected channels 
    bc = reshape(repmat(EEG.artifacts.BC,[1 nSm 1]),[nEl nSm*nEp]);
    bc = find(all(bc,2));
    chcol = repmat({'k'},[1 nEl]);
    for i=1:length(bc)
        chcol{bc(i)} = 'r';
    end
else
    chcol = {'k'};
end


eegplot(EEG.data,'srate',EEG.srate,'winlength',winlength,'color',chcol,'winrej',winrej);



end

