function writePDIPSamples(n, matDir)
    outputDir = "samples/samplesPDIP_N"+n+".txt";
    disp(['Guardando: '+ outputDir]);
    load(matDir);

    I = size(Mx,1);
    nSamples = size(u_tilde,3);

    fileID = fopen(outputDir,'w');
    fprintf(fileID,"N:\t%d\tI:\t%d\tnSamples:\t%d\n", n, I, nSamples);
    fprintf(fileID,'%.30f\t',reshape(H',1,[]));
    fprintf(fileID,'\n');    
    fprintf(fileID,'%.30f\t',reshape(Mx',1,[]));
    fprintf(fileID,'\n'); 

    for sample = 1:nSamples
        fprintf(fileID,'%.30f\t',reshape(h(:,:,sample)',1,[]));
        fprintf(fileID,'\n');
        fprintf(fileID,'%.30f\t',reshape(cx(:,:,sample)',1,[]));
        fprintf(fileID,'\n');
        fprintf(fileID,'%.30f\t',reshape(u_tilde(:,:,sample),1,[]));
        fprintf(fileID,'\n');
    end
    fclose(fileID);
end


