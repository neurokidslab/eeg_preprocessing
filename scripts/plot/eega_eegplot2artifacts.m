function EEG = eega_eegplot2artifacts(EEG, rej)

if ndims(EEG.data)==3
    error('The function only works for continuos data')
end

if ~isfield(EEG, 'artifacts')
    nEl = size(EEG.data,1);
    nS = size(EEG.data,2);
    nEp = size(EEG.data,3);

    EEG.artifacts.algorithm.parameters = [];
    EEG.artifacts.algorithm.stepname = [];
    EEG.artifacts.algorithm.rejxstep = [];
    EEG.artifacts.BCT = false(nEl,nS,nEp);
    EEG.artifacts.BC = false(nEl,1,nEp);
    EEG.artifacts.BCmanual = [];
    EEG.artifacts.BT = false(1,nS,nEp);
    EEG.artifacts.BE = false(1,1,nEp);
    EEG.artifacts.BS = false(1);
end

btall = false(size(EEG.artifacts.BT));
for i=1:size(rej,1)
    btall(1,int32(rej(i,1)):int32(rej(i,2))) = 1;
end
btnew = btall & ~EEG.artifacts.BT;
btremoved = ~btall & EEG.artifacts.BT;

EEG.artifacts.BT = btall;
EEG.artifacts.BCT(:,btnew) = 1;
EEG.artifacts.BCT(:,btremoved) = 0;
EEG.artifacts.BTmanual = btnew;

end