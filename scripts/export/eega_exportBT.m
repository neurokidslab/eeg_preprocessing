function eega_exportBT(EEG, filename)

bt = reshape(EEG.artifacts.BT,[1 size(EEG.artifacts.BT,2)*size(EEG.artifacts.BT,3)]);

%get bad segment
badi = find( [ bt(1) diff(bt)==1 ] )';
badf = find( [ diff(bt)==-1 bt(end) ] )';
duration = badf-badi;

%convert to seconds
badi = badi/EEG.srate;
badf = badf/EEG.srate;
duration = duration/EEG.srate;

% description
description = repmat({'bad_apice'},[length(badi) 1]);

% save in a text file
% filename = fullfile(outputfolder,[filename '_BT.txt']);
fid=fopen(filename,'w');
fprintf(fid,'onset, duration, description\n');
for i=1:length(badi)
    fprintf(fid,'%f,%f,%s\n',badi(i),duration(i),description{i});
end
fclose(fid);

end