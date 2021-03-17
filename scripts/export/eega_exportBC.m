function EEG = eega_exportBC(EEG,outputfolder)

if nargin<2
    outputfolder = pwd;
end

% badEpoch = EEG.artifacts.BE;
badCh = any(EEG.artifacts.BC,3);
% badCh = badCh(:,:,~badEpoch);


badCh = any(badCh,3);
isbadCh = find(badCh);

[~,filename,~] = fileparts(EEG.filename);
filename = fullfile(outputfolder,[filename '_BC.txt']);
fid=fopen(filename,'w');
for i=1:length(isbadCh)
    fprintf(fid,'%d\n',isbadCh(i));
end
fclose(fid);
end

