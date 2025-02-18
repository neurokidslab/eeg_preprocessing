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
        for k=1:length(sbjfile)
            if ~strcmp(sbjfile(k).name(1),'.')
                sbjfilek = fullfile(sbjfile(k).folder,sbjfile(k).name);
                fprintf('--> Adding file: %s \n',sbjfilek)
                SBJs = cat(1,SBJs,{sbjfilek});
            end
        end
    else
        fprintf('--> NO file was find in folder %s\n',sbjfolder)
    end
end


end