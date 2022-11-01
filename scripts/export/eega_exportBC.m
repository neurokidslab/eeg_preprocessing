function EEG = eega_exportBC(EEG,filename)

fid=fopen(filename,'w');
fprintf(fid,'# APICE bad channels\n');
fprintf(fid,'# bad channels per epoch (in raws)\n');

for ep=1:size(EEG.artifacts.BC,3)
    badCh = find(EEG.artifacts.BC(:,1,ep));
    badCh = find(badCh);

    % save in a text file
    for i=1:length(badCh)
        fprintf(fid,'%s,',EEG.chanlocs(badCh(i)).labels);
    end
    fprintf(fid,'\n');

end
fclose(fid);
end

