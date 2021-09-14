% this function uses eegplotto plot the rejection matrixes
% By default it plots BT and BC(all recording)
% It can also plot BTC (bad times x channles)
% this can be defined by setting the optional inputs 
% - plotbct (defaoult 0)
% - plotbt (defaoult 1)
% - plotbc (defaoult 1)

function eega_plot_rejection(EEG, plotbct, plotbt, plotbc, winlength)

if nargin<5 || isempty(winlength)
    if size(EEG.data,3)==1
        winlength = min(100, size(EEG.data,2)/EEG.srate);
    else
        winlength = 10;
    end
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

[nEl, nSm, nEp] = size(EEG.data);
if nEp>1 && plotbt
    error('This function does not support plotting BT for epoched data')
end

% create the matrix holding what to mark as rejected
winrej = [];
colrejbct = [1 1 1];
if plotbct  % plot in red the rejected data (white background)
    bct = reshape(EEG.artifacts.BCT,[nEl nSm*nEp]);
    for el=1:nEl
        btel = bct(el,:);
        btel_ini = find([btel(1) diff(btel,1,2)==1])';
        btel_fin = find([diff(btel,1,2)==-1 btel(end)])';
        ch = false(length(btel_ini), nEl);
        ch(:,el) = 1;
        winrejel = cat(2,btel_ini,btel_fin, repmat(colrejbct,[length(btel_ini) 1]), ch);
        winrej = cat(1,winrej,winrejel);
    end
end

colrejbt = [0.9 0.60 0.25];
if plotbt % plot with orange background the rejected times 
    bt = reshape(EEG.artifacts.BT,[1 nSm*nEp]);
    bt_ini = find([bt(1) diff(bt,1,2)==1])';
    bt_fin = find([diff(bt,1,2)==-1 bt(end)])';
%     if nEp>1
%         for ep=1:nEp
%             ini_ep = (ep-1)*nSm+1;
%             fin_ep = ep*nSm;
%             cut_idx = ( (bt_ini>ini_ep) & (bt_ini<fin_ep) );
%             if any(cut_idx)
%                 bt_ini = cat(1, bt_ini, fin_ep+1);
%                 bt_fin = cat(1, bt_fin, fin_ep);
%             end
%             bt_ini = sort(bt_ini);
%             bt_fin = sort(bt_fin);
%         end
%     end

    winrejel = cat(2,bt_ini,bt_fin, repmat(colrejbt,[length(bt_ini) 1]), true(length(bt_ini), nEl));
    winrej = cat(1,winrej,winrejel);
end

if plotbc % plot in red the rejected channels 
    % channels rejected during al the recording 
    bc = reshape(repmat(EEG.artifacts.BC,[1 nSm 1]),[nEl nSm*nEp]);
    bc = find(all(bc,2));
    chcol = repmat({'k'},[1 nEl]);
    for i=1:length(bc)
        chcol{bc(i)} = 'r';
    end
    % channele rejected per epoch
    if nEp>1
        for ep=1:nEp
            badch = find(EEG.artifacts.BC(:,1,ep));
            for i=1:length(badch)
                ch = false(1, nEl);
                ch(badch(i)) = 1;
                bteli_ini = (ep-1)*nSm+1;
                bteli_fin = ep*nSm;
                winreji = [bteli_ini bteli_fin 1 1 1 ch];
                winrej = cat(1,winrej,winreji);
            end
        end
    end
else
    chcol = {'k'};
end


eegplot(EEG.data,'srate',EEG.srate,'winlength',winlength,'color',chcol,'winrej',winrej);



end

