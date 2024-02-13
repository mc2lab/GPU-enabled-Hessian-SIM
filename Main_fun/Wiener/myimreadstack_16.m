function img_read=myimreadstack_16(filename,nk,zhens,testreadx,testready)
    warning off;
	t = Tiff(filename, 'r');
    t.setDirectory(nk);
    img_read=zeros(testreadx,testready,zhens,'uint16');
    for k = 1:zhens-1
		img_read(:,:,k)=t.read();
        t.nextDirectory();
    end
    img_read(:,:,k+1)=t.read();
	t.close();
    warning on;
end