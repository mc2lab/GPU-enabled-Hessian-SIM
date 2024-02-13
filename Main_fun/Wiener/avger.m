function avger(pathname,filename)
    starframe=1;
    avger_num=9;
    wenjianming=(fullfile(pathname, filename));
    info = imfinfo(wenjianming);
    zstack_all = numel(info);
    sizex=info(1).Height;
    sizey=info(1).Width;
    avger_result=zeros(sizex,sizey,ceil((zstack_all-(starframe-1))./avger_num),'uint16');
    clear info;
    fd_all=myimreadstack_16(wenjianming,1,zstack_all,sizex,sizey);
    for iiw=1:avger_num:(zstack_all-(avger_num-1)-(starframe-1))
        avger_tmp=sum(fd_all(:,:,iiw:iiw+(avger_num-1)),3,'native')./avger_num;
        avger_result(:,:,ceil((iiw+(avger_num-1))./(avger_num)))=avger_tmp;
    end
        imwritestack_16(avger_result,[ pathname  '\Pseudo-TIRF\TIRF_' filename(1:end-4) '.tif']);
end