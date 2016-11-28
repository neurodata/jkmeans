m = ones(28*28, 0);
for i = 0:9
    infile = strcat('data',num2str(i))
    fid = fopen(infile);
    [ t1, N ] = fread(fid, [28 * 28 1000], 'uchar');
    fclose(fid);
    m = [m t1];
end

ofile = strcat('mnist.csv');
csvwrite(ofile, m);