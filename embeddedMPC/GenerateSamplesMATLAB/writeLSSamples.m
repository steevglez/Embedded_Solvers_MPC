function writeLSSamples(n, matDir)

    outputDir = "samples/samplesLS_N"+n+".txt";
    disp(['Guardando: '+ outputDir]);
    load(matDir);

    nSamples = size(Ak,3);
    fileID = fopen(outputDir,'w');
    fprintf(fileID,"N:\t%d\tI:\t%d\tnSamples:\t%d\n", n, nSamples);

    for sample = 1:nSamples
        fprintf(fileID,'%.30f\t',reshape(Ak(:,:,sample)',1,[]));
        fprintf(fileID,'\n');
        fprintf(fileID,'%.30f\t',reshape(bk(:,:,sample)',1,[]));
        fprintf(fileID,'\n');
        fprintf(fileID,'%.30f\t',reshape(zk(:,:,sample),1,[]));
        fprintf(fileID,'\n');
    end
    fclose(fileID);
end

