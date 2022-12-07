function SBJs = eega_getfilesinfolders(Path2Data, FilesIn)

% Generate a list of files to analyze 
sbjs2run = dir(Path2Data);
sbjs2run = sbjs2run(~ismember({sbjs2run.name},{'.', '..'}));
sbjs2run = sbjs2run([sbjs2run.isdir]);
SBJs = {};
for sbj=1:length(sbjs2run)
    
    sbjfolder = fullfile(sbjs2run(sbj).folder,sbjs2run(sbj).name);
    sbjfile = dir(fullfile(sbjfolder,FilesIn));
    if ~isempty(sbjfile)
        sbjfile = fullfile(sbjfile(sbj).folder,sbjfile(sbj).name);
        fprintf('--> Adding file: %s \n',sbjfile)
        SBJs = cat(1,SBJs,{sbjfile});
    else
        fprintf('--> NO file was find in folder %s\n',sbjfolder)
    end
end


end