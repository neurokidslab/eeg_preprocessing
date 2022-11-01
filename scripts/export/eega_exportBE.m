function EEG = eega_exportBE(EEG,filename)

badEpoch = find(EEG.artifacts.BE);

fid=fopen(filename,'w');
fprintf(fid,'# APICE bad epochs\n');
for i=1:length(badEpoch)
    fprintf(fid,'%d\n',badEpoch(i));
end
fclose(fid);
end

