function EEG = eega_exportBE(EEG,outputfolder)

if nargin<2
    outputfolder = pwd;
end

badEpoch = find(EEG.artifacts.BE);

[~,filename,~] = fileparts(EEG.filename);
filename = fullfile(outputfolder,[filename '_BE.txt']);
fid=fopen(filename,'w');
for i=1:length(badEpoch)
    fprintf(fid,'%d\n',badEpoch(i));
end
fclose(fid);
end

