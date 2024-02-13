    %all
    if ((exist('flag','var')))
        clearvars -except wenjianming pathname filename psf n weilac spjg flag bgflag bgname starframe beishu_an beishu_re filename_all filei wavelengh pg fanwei regul deconv_otfname zstack testreadx testready savec sigma mu nphases nangles Progressbar
        all_name=wenjianming;
        clear flag
    else
        error('Please choose the main function to Reconstruction')
    end
    phase_matrix= [1 1 1;exp(1i*regul*(spjg(1)/sum(spjg))) 1  exp(-1i*regul*(spjg(1)/sum(spjg)));exp(1i*regul*((spjg(1)+spjg(2))/sum(spjg))) 1  exp(-1i*regul*((spjg(1)+spjg(2))/sum(spjg)));];%488

    
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Pre-process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    n_512=max([512,testreadx,testready]);
    fc_ang = ceil(120*(n_512/512));
    if wavelengh==647
        fc_con = ceil(80*(n_512/512));
    else
        fc_con = ceil(105*(n_512/512));
    end
    jdjd=0.02;
    jdc=0.02;
    nphases = 3; 
    load  ([pathname  filename(1:end-4) '_zuobiaox.mat']);
    load  ([pathname  filename(1:end-4) '_zuobiaoy.mat']);
    zuobiaox_512=zuobiaox*(n_512/floor((zuobiaox(1,1)+zuobiaoy(1,1))/2));
    zuobiaoy_512=zuobiaoy*(n_512/floor((zuobiaox(1,1)+zuobiaoy(1,1))/2));
    H_ang = imreadstack(deconv_otfname);
    H_ang = imresize(H_ang,[n_512,n_512] ,'bilinear') ;
    H_con = H_ang;
    [k_x, k_y]=meshgrid(-(n_512)/2:(n_512)/2-1, -(n_512)/2:(n_512)/2-1);
    k_r = sqrt(k_x.^2+k_y.^2);
    indi_ang =  k_r > fc_ang ;
    indi_con =  k_r > fc_con ;
    H_ang(indi_ang)=0;
    H_con(indi_con)=0;
    H_ang=abs(H_ang);
    H_con=abs(H_con);
    H1_ang=H_ang;
    H1_con=H_con;
    H1_ang(H1_ang~=0)=1;
    H1_con(H1_con~=0)=1;
    H9_ang=repmat(H1_ang,[1 1 (nphases*nangles)]);
    H9_con=repmat(H1_con,[1 1 (nphases*nangles)]);
    if (starframe+(nphases*nangles)*beishu_an-1<=zstack)
        fd=myimreadstack_TIRF(all_name,starframe,(nphases*nangles)*beishu_an,testreadx,testready);
    else
        close(Progressbar);
        warndlg('The number of averaged frames should be smaller','warn','modal'); 
        return;
    end
    for beishui=1:1:(nphases*nangles)
        fdd_an(:,:,beishui)=sum(fd(:,:,[beishui:(nphases*nangles):(nphases*nangles)*beishu_an]),3)./beishu_an;
    end
    fd=fdd_an;
    fd512=zeros(n_512,n_512);
	[sx,sy,~]=size(fd);
    K_h = [size(fd,1),size(fd,2)];
    N_h = [n_512,n_512];
    L_h = ceil((N_h-K_h) / 2);
    v_h = colonvec(L_h+1, L_h+K_h);
    hw=zeros(N_h);
    if bgflag==1
        bg=imreadstack(bgname);
%         bgjie488=bg(1:sx,1:sy); 
        bgjie488=bg((end/2)+1-(sx)/2:(end/2)+(sx)/2,(end/2)+1-(sy)/2:(end/2)+(sy)/2); 
    elseif bgflag==0
        bgjie488=zeros(sx,sy); 
    end
    for ii=1:1:(nphases*nangles)
        fd(:,:,ii)=fd(:,:,ii)-bgjie488;
        hw(v_h{:})=fd(:,:,ii);
        fd512(:,:,ii)=hw;
    end
    DIbars = fftshift(fft2(ifftshift(fd512)));
    DIbars_ang = H9_ang .* DIbars;
    DIbars_con = H9_con .* DIbars;
    inv_phase_matrix = inv(phase_matrix);
    
    c6=zeros(1,1,3*(nphases-1));
    angle6=zeros(1,1,3*(nphases-1));
    replc6_con=zeros(2*n_512,2*n_512,(nangles)*(nphases-1));
    replc6_ang=zeros(2*n_512,2*n_512,(nangles)*(nphases-1));
    cm_con=zeros(2*n_512,2*n_512,(nangles)*(nphases-1));
    cm_ang=zeros(2*n_512,2*n_512,(nangles)*(nphases-1));
    R2_ang=zeros(2*n_512,2*n_512,(nangles)*(nphases-1));
    R2_zhongxin_ang=zeros(2*n_512,2*n_512,(nangles)*(nphases-1));
    sp = zeros(n_512,n_512,nphases*nangles);
    temp_separated = zeros(n_512,n_512,nphases); 
    x=0:(2*n_512-1);
    y=(0:(2*n_512-1))';
    xx2=repmat(x,2*n_512,1);
    yy2=repmat(y,1,2*n_512);
    
    
    K_h = size(H1_ang);
    N_h = 2*K_h;
    L_h = ceil((N_h-K_h) / 2);
    v_h = colonvec(L_h+1, L_h+K_h);
    hw=zeros(N_h);
    hw(v_h{:})=H_ang;
    H_ang=hw;
    hw(v_h{:})=H_con;
    H_con=hw;
    hw(v_h{:})=H1_ang;
    H1_ang=hw;
    hw(v_h{:})=H1_con;
    H1_con=hw;
    for itheta=1:nangles
         for j = 1:nphases
                for k = 1:nphases
                    temp_separated(:,:,k) = inv_phase_matrix(j,k).*DIbars(:,:,(itheta-1)*3+k);
                    sp(:,:,(itheta-1)*3+j) = sp(:,:,(itheta-1)*3+j)+temp_separated(:,:,k);
                end

         end
    end
    for spi=1:1:(nangles)*(nphases-1)
        spzhongxin = sp(:,:,ceil(spi/2)*3-1);
        spzhongyiwei = sp(:,:,ceil(spi/2)+2*floor(spi/2));
        hw(v_h{:})=spzhongyiwei;
        spzhongyiwei=hw;
        hw(v_h{:})=spzhongxin;
        spzhongxin=hw;
        kytest = 2*pi*(zuobiaox_512(spi+1+floor((spi-1)/2))-n_512)/(2*n_512);
        kxtest = 2*pi*(zuobiaoy_512(spi+1+floor((spi-1)/2))-n_512)/(2*n_512);
        Irtest(:,:)=exp(1i*(kxtest*xx2+kytest*yy2));
        replcHtest_ang(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1_ang))).*Irtest(:,:))));
        replcHtest_con(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1_con))).*Irtest(:,:))));
        replcHtest_ang(abs(replcHtest_ang)>0.9)=1;
        replcHtest_ang(abs(replcHtest_ang)~=1)=0;
        replcHtest_con(abs(replcHtest_con)>0.9)=1;
        replcHtest_con(abs(replcHtest_con)~=1)=0;
        replch_ang(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H_ang))).*Irtest(:,:))));
        replch_con(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H_con))).*Irtest(:,:))));
        replch_ang = replch_ang.*replcHtest_ang;
        replch_con = replch_con.*replcHtest_con;
        replctest(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(spzhongyiwei))).*Irtest(:,:))));
        youhuatest_ang=replctest.*replcHtest_ang.*H_ang;
        youhuatest_con=replctest.*replcHtest_con.*H_con;
        chongdie_ang=spzhongxin.*replch_ang;
        chongdie_con=spzhongxin.*replch_con;
        cm_con(:,:,spi)=youhuatest_con./(chongdie_con+eps);
        replc6_con(:,:,spi)=replch_con.*H_con;
        cm_ang(:,:,spi)=youhuatest_ang./(chongdie_ang+eps);
        replc6_ang(:,:,spi)=replch_ang.*H_ang;
        R2_ang(:,:,spi)=youhuatest_ang;
        R2_zhongxin_ang(:,:,spi)=chongdie_ang;
    end
    angcm=angle(cm_ang);
    abscm=abs(cm_con);
%     imwritestack((angcm),'angcm.tif');
%     imwritestack((abscm),'abscm.tif');
    f=-pi:jdjd:pi;
    ff=0:jdc:0.6;
    for ii=1:1:(nangles)*(nphases-1)
        a_ang=replc6_ang(:,:,ii);
        a_con=replc6_con(:,:,ii);
        b_ang=a_ang(:);
        b_con=a_con(:);
        b_ang(b_ang~=0)=1;
        b_con(b_con~=0)=1;
        c=angcm(:,:,ii);
        cc=abscm(:,:,ii);
        d=c(:);
        dd=cc(:);
        b_ang(:,2)=d;
        b_con(:,2)=dd;
        b_ang((b_ang(:,1)==0),:)=[] ;
        b_con((b_con(:,1)==0),:)=[] ;
        e=b_ang(:,2);
        ee=b_con(:,2);
        e=e';
        ee=ee';
        g=histc(e,f);
        gg=histc(ee,ff);
%         figure;plot(g);
%         hold on
        imf=emd(g);
        imf2=emd(gg);
		if max(g)<50
            g=sum(imf(5:end,:),1);
        else
            g=sum(imf(4:end,:),1);
        end
        gg=sum(imf2(1:end,:),1);
        p=g;
        pp=gg;
%         plot(p,'r');
%         hold off;
        h=find(g==max(g));
        hh=(find(gg==max(gg)));
        angle6(:,:,ii)=-pi+h*jdjd;
        c6(:,:,ii)=jdc*mean(hh);
        %c6(:,:,ii)=0.7;
        regress_R2;
    end
    clearvars -except c6 all_name zstack angle6 fenshu n dizhi  beishu_re weilac yupan sx sy starframe  deconv_otfname bgjie488 spjg pathname filename fd filename_all filei pg fanwei wavelengh regul testreadx testready savec sigma mu nphases nangles Progressbar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Pre-process
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Wiener SIM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    waitbar(0, Progressbar, 'Wiener Reconstruction.');
    load  ([pathname  filename(1:end-4) '_zuobiaox.mat']);
    load  ([pathname  filename(1:end-4) '_zuobiaoy.mat']);
    if savec==1
        save([pathname  filename(1:end-4) '_samec_angle.mat'],'c6','angle6');
    end
    image_size = max([sx,sy]);
    if image_size<=256
        otfx=256;
    elseif image_size>256 && image_size<=512
        otfx=512;
    elseif image_size>512
        otfx=image_size;
%         error('Image size could not larger than 512*512')
    end
    zuobiaox=zuobiaox*(otfx/n);
    zuobiaoy=zuobiaoy*(otfx/n);
    n=otfx;
    num_images=zstack-starframe+1;
    rep_images=10;
    xishu=[1,((1/c6(:,:,1))+(1/c6(:,:,2)))/2,((1/c6(:,:,1))+(1/c6(:,:,2)))/2,1,((1/c6(:,:,3))+(1/c6(:,:,4)))/2,((1/c6(:,:,3))+(1/c6(:,:,4)))/2,1,((1/c6(:,:,5))+(1/c6(:,:,6)))/2,((1/c6(:,:,5))+(1/c6(:,:,6)))/2];
    modulation_depth=xishu
    fc = ceil(200*(n/512)); % maximum spatial frequency
    plong= floor(sum(((zuobiaox-zuobiaox(1,1)).^2+(zuobiaoy-zuobiaoy(1,1)).^2).^0.5)./(2*nangles));
    kekao=ceil(135*(n/512)); 
    deph1= sign(angle6(:,:,1))*(abs(angle6(:,:,1))+abs(angle6(:,:,2)))*0.5;
    deph2= sign(angle6(:,:,3))*(abs(angle6(:,:,3))+abs(angle6(:,:,4)))*0.5;
    deph3= sign(angle6(:,:,5))*(abs(angle6(:,:,5))+abs(angle6(:,:,6)))*0.5;
    fd512=zeros(n,n);
    tirff=zeros(2*sx,2*sy,min(floor(ceil(num_images)/((nphases*nangles)*beishu_re)),rep_images),'single');
    zhen=starframe;
    DIbar = zeros(n,n,(nphases*nangles));
    retirff = zeros(2*n,2*n,(nphases*nangles));
    spsim = zeros((2*n),(2*n),nphases*nangles);
    replcHtest = zeros((2*n),(2*n),(nphases*nangles));
    replch = zeros((2*n),(2*n),(nphases*nangles));
    tmprc1 = zeros((2*n),(2*n),(nphases*nangles));
    Irtest = zeros((2*n),(2*n),(nphases*nangles));
    x=0:(2*n-1);
    y=(0:(2*n-1))';
    xx2=repmat(x,2*n,1);
    yy2=repmat(y,1,2*n);
    
    [k_x, k_y]=meshgrid(-(n)/2:(n)/2-1, -(n)/2:(n)/2-1);
    k_r = sqrt(k_x.^2+k_y.^2);
    indi =  k_r > fc ;
    jiequ=ones(n,n);
    jiequ(indi) = 0;
    jiequ9=repmat(jiequ,[1,1,(nphases*nangles)]);
    otf = imreadstack(deconv_otfname);
    if size(otf,1)~=n || size(otf,2)~=n
        psf = imresize(otf,[n,n] ,'bilinear') ;
    else
        psf = otf;
    end
    H = psf;
    H=H.*jiequ;
    H=H./max(H(:));
    H1=H;
    H1(H1~=0)=1;
    K_h = size(psf);
    N_h = 2*K_h;
    L_h = ceil((N_h-K_h) / 2);
    v_h = colonvec(L_h+1, L_h+K_h);
    hw=zeros(N_h);
    hw(v_h{:})=H;
    Hk=hw;
    hw(v_h{:})=H1;
    H1=hw;
    clear hw N_h K_h L_h v_h
	for ii=1:1:(nphases*nangles)
		kytest = 2*pi*(zuobiaox(ii,:)-n)/(2*n);
		kxtest = 2*pi*(zuobiaoy(ii,:)-n)/(2*n);
		Irtest(:,:,ii)=exp(1i*(kxtest*xx2+kytest*yy2));
		replcHtest(:,:,ii) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1))).*Irtest(:,:,ii))));
		replch(:,:,ii) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(Hk))).*Irtest(:,:,ii))));
	end
	replcHtest(abs(replcHtest)>0.9)=1;
	replcHtest(abs(replcHtest)~=1)=0;
	replch= replch.*replcHtest; 
	re=replch;
	absre=abs(re);
	absre(absre~=0)=1;
	re=re.*absre;
	re(re==0)=10000000000000000000;
	re=re./abs(re);
	reh=abs(replch);
	hs = sum(abs(reh(:,:,:)).^2,3);
	reh=repmat(reh,[1 1 2]);
	[k_x, k_y]=meshgrid(-(2*n)/2:(2*n)/2-1, -(2*n)/2:(2*n)/2-1);
	k_r = sqrt(k_x.^2+k_y.^2);
	k_max = plong+fc;
	bhs = cos(pi*k_r/(2*k_max));
	indi = find( k_r > k_max );
	bhs(indi) = 0;
	clear H1 H2 Hk absre;
	
	phi1tmp=linspace(0+deph1,regul+deph1,sum(spjg)+1);
    phi1(1)=phi1tmp(1);
    phi1(2)=phi1tmp(1+spjg(1));
    phi1(3)=phi1tmp(1+spjg(1)+spjg(2));

	
	phi2tmp=linspace(0+deph2,regul+deph2,sum(spjg)+1);
    phi2(1)=phi2tmp(1);
    phi2(2)=phi2tmp(1+spjg(1));
    phi2(3)=phi2tmp(1+spjg(1)+spjg(2));



	phi3tmp=linspace(0+deph3,regul+deph3,sum(spjg)+1);
    phi3(1)=phi3tmp(1);
    phi3(2)=phi3tmp(1+spjg(1));
    phi3(3)=phi3tmp(1+spjg(1)+spjg(2));

	phase_matrix1 = [1 exp(1i*phi1(1)) exp(-1i*phi1(1));1 exp(1i*phi1(2)) exp(-1i*phi1(2));1 exp(1i*phi1(3)) exp(-1i*phi1(3));];
	phase_matrix2 = [1 exp(1i*phi2(1)) exp(-1i*phi2(1));1 exp(1i*phi2(2)) exp(-1i*phi2(2));1 exp(1i*phi2(3)) exp(-1i*phi2(3));];
	phase_matrix3 = [1 exp(1i*phi3(1)) exp(-1i*phi3(1));1 exp(1i*phi3(2)) exp(-1i*phi3(2));1 exp(1i*phi3(3)) exp(-1i*phi3(3));];
	inv_phase_matrix1 = inv(phase_matrix1);
	inv_phase_matrix2 = inv(phase_matrix2);
	inv_phase_matrix3 = inv(phase_matrix3);
	inv_phase_matrix(:,:,1) = inv_phase_matrix1;
	inv_phase_matrix(:,:,2) = inv_phase_matrix2;
	inv_phase_matrix(:,:,3) = inv_phase_matrix3;
	
	padsize=0;
    x = 1:(sy+2*padsize);
    y = (1:(sx+2*padsize))';
    sig=0.25;
    mask = repmat(sigmoid(sig*(x-padsize)) - sigmoid(sig*(x-sy-padsize-1)), sx+2*padsize, 1) .* repmat(sigmoid(sig*(y-padsize)) - sigmoid(sig*(y-sx-padsize-1)), 1, sy+2*padsize);
    mask9=repmat(mask.^3,[1 1 (nphases*nangles)]);
    K_h2 = [n,n];
    N_h2 = 2*K_h2;
    L_h2 = ceil((N_h2-K_h2) / 2);
    v_h2 = colonvec(L_h2+1, L_h2+K_h2);
    hw2=zeros(N_h2);
while( zhen<=(num_images-(nphases*nangles)*beishu_re+starframe))
    if (zhen+(nphases*nangles)*beishu_re-1<=zstack)
        D=myimreadstack_TIRF(all_name,zhen,(nphases*nangles)*beishu_re,testreadx,testready); 
    else
        error('The parameter beishu_re should be smaller');
    end
    D=D-repmat(bgjie488,[1,1,(nphases*nangles)*beishu_re]);
    for beishui=1:1:(nphases*nangles)
        fdd_an(:,:,beishui)=sum(D(:,:,[beishui:(nphases*nangles):(nphases*nangles)*beishu_re]),3)./beishu_re;
    end
    D=fdd_an;
%     weights =mean(mean(double(D))); 
    %[num,val]=sort(weights);
    %paixu=(val(:,:,7:9));
    %panduan=ceil((paixu-0.1)./3);
    %if (any(diff(panduan)))||any(panduan~=yupan)||((num(9)-num(8))>(num(7)-num(6)))||((num(8)-num(7))>(num(7)-num(6)))||((num(9)-num(7))>(num(7)-num(6)))
    %    zhen=zhen+1;
    %    continue;
    %end
%     weights = mean(weights)./weights;
    D=D.*mask9;
    
    K_h = [size(D,1),size(D,2)];
    N_h = [n,n];
    L_h = ceil((N_h-K_h) / 2);
    v_h = colonvec(L_h+1, L_h+K_h);
    hw=zeros(N_h);
    for I=1:(nphases*nangles)
%         D(:,:,I) = D(:,:,I) .* weights(1,1,I); 
        hw(v_h{:})=D(:,:,I);
        fd512(:,:,I)=hw;
        DIbar(:,:,I) = fftshift(fft2(ifftshift(fd512(:,:,I))));
    end

    clear D fd512;
        sp = zeros(n,n,nphases*nangles);
        for itheta=1:nangles
            for j = 1:nphases
                temp_separated = zeros(n,n,nphases);
                for k = 1:nphases
                    temp_separated(:,:,k) = inv_phase_matrix(j,k,itheta).*DIbar(:,:,(itheta-1)*3+k);
                    sp(:,:,(itheta-1)*3+j) = sp(:,:,(itheta-1)*3+j)+temp_separated(:,:,k);
                end
            end
        end
        sp=sp.*jiequ9;
        for ii=1:1:(nphases*nangles)
            hw2(v_h2{:})=sp(:,:,ii);
            spsim(:,:,ii)=hw2;
            retirff(:,:,ii)= fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(spsim(:,:,ii)))).*Irtest(:,:,ii))));
        end

        spttH=retirff.*replcHtest;
        retirff=spttH./(re);
        
        for t = 1:nphases*nangles
            tmprc1(:,:,t) =   (xishu(t).*retirff(:,:,t).*(conj(reh(:,:,t))))./ ( hs + .005*length(itheta)*(weilac)^2);  
        end
		dr = sum(tmprc1,3);
        drr = dr.*bhs;
        fimage = fftshift(ifft2(ifftshift(drr)));
        tirff(:,:,1+mod(-1+round(ceil(zhen-starframe+(nphases*nangles)*beishu_re)/((nphases*nangles)*beishu_re)),rep_images)) = fimage((end/2)+1-(sx):(end/2)+(sx),(end/2)+1-(sy):(end/2)+(sy));
        if 1+mod(-1+round(ceil(zhen-starframe+(nphases*nangles)*beishu_re)/((nphases*nangles)*beishu_re)),rep_images)==rep_images
            if round(ceil(zhen-starframe+(nphases*nangles)*beishu_re)/((nphases*nangles)*beishu_re))==rep_images
                tirff(tirff<0)=0;
                imwritestack(tirff,[pathname   '\SIM-Wiener\re-' filename]);  
            else
                tirff(tirff<0)=0;
                imwritestack_appe(tirff,[pathname   '\SIM-Wiener\re-' filename]);  
            end
        end
        zhen=zhen+(nphases*nangles)*beishu_re;
        waitbar(zhen/(num_images-(nphases*nangles)*beishu_re+starframe), Progressbar, 'Wiener Reconstruction.');
end
if 1+mod(-1+round(ceil(zhen-starframe+(nphases*nangles)*beishu_re)/((nphases*nangles)*beishu_re)),rep_images) ~= 1
    if round(ceil(zhen-starframe+(nphases*nangles)*beishu_re)/((nphases*nangles)*beishu_re)) <= rep_images
        tirff(tirff<0)=0;
        imwritestack(tirff(:,:,1:mod(-1+round(ceil(zhen-starframe+(nphases*nangles)*beishu_re)/((nphases*nangles)*beishu_re)),rep_images)),[pathname   '\SIM-Wiener\re-' filename]);  
    else
        tirff(tirff<0)=0;
        imwritestack_appe(tirff(:,:,1:mod(-1+round(ceil(zhen-starframe+(nphases*nangles)*beishu_re)/((nphases*nangles)*beishu_re)),rep_images)),[pathname   '\SIM-Wiener\re-' filename]);  
    end
end
disp('Wiener reconstruction Successfully');
waitbar(0, Progressbar, 'Hessian reconstruction');
Bregman_Hessian_LowRam_2;
waitbar(0, Progressbar, 'TV reconstruction');
Bregman_TV_denoise;
close(Progressbar);
Running_average;
helpdlg('All Done');