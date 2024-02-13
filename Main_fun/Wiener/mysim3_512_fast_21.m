clearvars -except pathname filename psf spjg weilac bgflag bgname starframe beishu_an beishu_re filename_all filei fanwei pg wavelengh regul deconv_otfname zstack testreadx testready savec notch_swith notch_para_1 notch_para_2 sigma mu nphases nangles Progressbar
wenjianming=(fullfile(pathname, filename));
if (starframe+(nphases*nangles)*beishu_an-1<=zstack)
    fd=myimreadstack_TIRF(wenjianming,starframe,(nphases*nangles)*beishu_an,testreadx,testready);
else
    close(Progressbar);
    warndlg('The number of averaged frames should be smaller','warn','modal'); 
    return;
end
%fd=fd(:,:,19:27);
psf = imresize(psf,[512,512] ,'bilinear') ;%add
H = psf;
[fdx,fdy,fdz]=size(fd);
fdd=zeros(fdx,fdy,(nphases*nangles));
n = max([fdx,fdy,512]);
fc = ceil(220*(n/512));
pg = ceil(pg*(n/512));
fanwei = ceil(fanwei*(n/512));
if n>512
     H = imresize(H,[n,n] ,'bilinear') ;
end
phase_matrix= [1 1 1;exp(1i*regul*(spjg(1)/sum(spjg))) 1  exp(-1i*regul*(spjg(1)/sum(spjg)));exp(1i*regul*((spjg(1)+spjg(2))/sum(spjg))) 1  exp(-1i*regul*((spjg(1)+spjg(2))/sum(spjg)));];%488
% padsize=0;
% xmask = 1:(fdy+2*padsize);
% ymask = (1:(fdx+2*padsize))';
% sigma=0.25;
% mask = repmat(sigmoid(sigma*(xmask-padsize)) - sigmoid(sigma*(xmask-fdy-padsize-1)), fdx+2*padsize, 1) .* repmat(sigmoid(sigma*(ymask-padsize)) - sigmoid(sigma*(ymask-fdx-padsize-1)), 1, fdy+2*padsize);
% mask9=repmat(mask.^3,[1 1 (nphases*nangles)]);

for beishui=1:1:(nphases*nangles)
    fdd(:,:,beishui)=sum(fd(:,:,[beishui:(nphases*nangles):(nphases*nangles)*beishu_an]),3)./beishu_an;
end
% fd=fdd.*mask9;
fd=fdd;
jdjd=0.02;
test=zeros(4,1);
zuobiaox(1:(nphases*nangles),1)=n;
zuobiaoy(1:(nphases*nangles),1)=n;
cm=zeros(2*n,2*n,(nangles)*(nphases-1));
replc6=zeros(2*n-1,2*n-1,(nangles)*(nphases-1));
cmabs=zeros(1,1,(nangles)*(nphases-1));
cmang=zeros(1,1,(nangles)*(nphases-1));
angle6=zeros(1,1,(nangles)*(nphases-1));
[k_x, k_y]=meshgrid(-(n)/2:(n)/2-1, -(n)/2:(n)/2-1);
k_r = sqrt(k_x.^2+k_y.^2);
indi =  k_r > fc ;
H(indi)=0;
H=abs(H);
inv_phase_matrix = inv(phase_matrix);
temp_separated = zeros(n,n,nphases);
sp = zeros(n,n,nphases*nangles);
spt = zeros(2*n,2*n,nphases*nangles);

fd512=zeros(n,n);
K_h = [size(fd,1),size(fd,2)];
N_h = [n,n];
L_h = ceil((N_h-K_h) / 2);
v_h = colonvec(L_h+1, L_h+K_h);
hw=zeros(N_h);
for ii=1:1:(nphases*nangles)
    hw(v_h{:})=fd(:,:,ii);
    fd512(:,:,ii)=hw;
end
DIbars = fftshift(fft2(ifftshift(fd512)));

H1=H;
H1(H1~=0)=1;
H2=H1;
H9=repmat(H1,[1 1 (nphases*nangles)]);
DIbars = H9 .* DIbars;

K_h = size(H2);
N_h = 2*K_h;
L_h = ceil((N_h-K_h) / 2);
v_h = colonvec(L_h+1, L_h+K_h);
hw=zeros(N_h);
hw(v_h{:})=H;
H=hw;
hw(v_h{:})=H1;
H1=hw;
for itheta=1:nangles
     for j = 1:nphases
            for k = 1:nphases
                temp_separated(:,:,k) = inv_phase_matrix(j,k).*DIbars(:,:,(itheta-1)*nphases+k);
                sp(:,:,(itheta-1)*nphases+j) = sp(:,:,(itheta-1)*nphases+j)+temp_separated(:,:,k);
            end

     end
end
sp=sp./(abs(sp)+eps);
for spi=1:2:(nangles)*(nphases-1)
 
 buchangx=1;
 buchangy=1;
 sp = sp .* H9; 
 spzhongxin = sp(:,:,ceil(spi/2)*nphases-1);
 spzhongyiwei = sp(:,:,ceil(spi/2)+2*floor(spi/2));
 Hzuan = H2(end:-1:1,end:-1:1);
 cishu=dft(H2,Hzuan);
 ci=cishu;
 ci(ci<0.9)=0;
 ci(ci~=0)=1;
 sp_tmp=conj(spzhongyiwei);
 sp_tmp=sp_tmp(end:-1:1,end:-1:1);
 jieguo=dft(spzhongxin,sp_tmp);
 jieguo=jieguo.*ci;
 lihe=abs(jieguo./(cishu+eps));
 [k_x, k_y]=meshgrid(-(2*n)/2+1:(2*n)/2-1, -(2*n)/2+1:(2*n)/2-1);
 k_r = sqrt(k_x.^2+k_y.^2);
 indi =  ((k_r < (pg-fanwei))|(k_r > (pg+fanwei))) ;
 lihe(indi)=0;
 [xx,yy]=find(lihe==max(max((lihe))));

 
 
 maxx=xx;
 maxy=yy;
 ky = 2*pi*(maxx-n)/(2*n);
 kx = 2*pi*(maxy-n)/(2*n);
 x=0:(2*n-1);
 y=(0:(2*n-1))';
 xx2=repmat(x,2*n,1);
 yy2=repmat(y,1,2*n);
 

hw(v_h{:})=spzhongyiwei;
spzhongyiwei=hw;
hw(v_h{:})=spzhongxin;
spzhongxin=hw;

 Irtest(:,:)=exp(1i*(kx*xx2+ky*yy2));
 replcHtest(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1))).*Irtest(:,:))));
 replcHtest(abs(replcHtest)>0.9)=1;
 replcHtest(abs(replcHtest)~=1)=0;
 
 replch(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H))).*Irtest(:,:))));
 replch= replch.*replcHtest;
 replctest(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(spzhongyiwei))).*Irtest(:,:))));
 youhuatest=replctest.*replcHtest.*H;
 hetest=(conj(youhuatest)).*(spzhongxin.*replch);
 hetest=abs(sum(hetest(:)));
 cishuhtest=H1.*replcHtest;
 cishuhtest=sum(cishuhtest(:));
 hetest=hetest./cishuhtest;
 he=hetest;
 


 maxx_tmp1=maxx-10^-5;
 maxx_tmp2=maxx+10^-5;
 maxy_tmp1=maxy-10^-5;
 maxy_tmp2=maxy+10^-5;
 ky_tmp1 = 2*pi*(maxx_tmp1-n)/(2*n);
 kx_tmp1 = 2*pi*(maxy_tmp1-n)/(2*n);
 ky_tmp2 = 2*pi*(maxx_tmp2-n)/(2*n);
 kx_tmp2 = 2*pi*(maxy_tmp2-n)/(2*n);
 for ii=1:1:4
 switch ii
     case 1
         kxtest=kx_tmp1;
         kytest=ky;
     case 2
         kxtest=kx_tmp2;
         kytest=ky;
     case 3
         kxtest=kx;
         kytest=ky_tmp1;
     case 4
         kxtest=kx;
         kytest=ky_tmp2;
 end
 Irtest(:,:)=exp(1i*(kxtest*xx2+kytest*yy2));
 replcHtest(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1))).*Irtest(:,:))));
 replcHtest(abs(replcHtest)>0.9)=1;
 replcHtest(abs(replcHtest)~=1)=0;
 
 replch(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H))).*Irtest(:,:))));
 replch= replch.*replcHtest;
 replctest(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(spzhongyiwei))).*Irtest(:,:))));
 youhuatest=replctest.*replcHtest.*H;
 hetest=(conj(youhuatest)).*(spzhongxin.*replch);
 hetest=abs(sum(hetest(:)));
 cishuhtest=H1.*replcHtest;
 cishuhtest=sum(cishuhtest(:));
 hetest=hetest./cishuhtest;
 test(ii)=hetest;
 end
 if((test(1)>test(2)))
     flag_maxy=-1;
 elseif((test(1)<test(2)))
     flag_maxy=+1;
 else
    close(Progressbar);
    warndlg('Can not estimate the pattern wave vector','warn','modal'); 
    return;
 end
 
 if((test(3)>test(4)))
     flag_maxx=-1;
 elseif((test(3)<test(4)))
     flag_maxx=+1;
 else
    close(Progressbar);
    warndlg('Can not estimate the pattern wave vector','warn','modal'); 
    return;
 end 
 maxx_tmp=maxx;
 maxy_tmp=maxy;
while ((buchangx>(10^-4))||(buchangy>(10^-4)))
%%%%%%%%%%%%%%%%%%%%%%%maxx%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
             maxx_tmp1=maxx-10^-5;
             maxx_tmp2=maxx+10^-5;
             ky_tmp1 = 2*pi*(maxx_tmp1-n)/(2*n);
             ky_tmp2 = 2*pi*(maxx_tmp2-n)/(2*n);
             for ii=3:1:4
                 switch ii
                     case 3
                         kxtest=2*pi*(maxy-n)/(2*n);
                         kytest=ky_tmp1;
                     case 4
                         kxtest=2*pi*(maxy-n)/(2*n);
                         kytest=ky_tmp2;
                 end
                 Irtest(:,:)=exp(1i*(kxtest*xx2+kytest*yy2));
                 replcHtest(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1))).*Irtest(:,:))));
                 replcHtest(abs(replcHtest)>0.9)=1;
                 replcHtest(abs(replcHtest)~=1)=0;

                 replch(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H))).*Irtest(:,:))));
                 replch= replch.*replcHtest;
                 replctest(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(spzhongyiwei))).*Irtest(:,:))));
                 youhuatest=replctest.*replcHtest.*H;
                 hetest=(conj(youhuatest)).*(spzhongxin.*replch);
                 hetest=abs(sum(hetest(:)));
                 cishuhtest=H1.*replcHtest;
                 cishuhtest=sum(cishuhtest(:));
                 hetest=hetest./cishuhtest;
                 test(ii)=hetest;
             end                 
             if((test(3)>test(4)))
                 flag_maxx=-1;
             elseif((test(3)<test(4)))
                 flag_maxx=+1;
             else
                 %error('Can not estimate the pattern wave vector');
                 flag_maxx=-1* flag_maxx;
             end    
     while(buchangx>(10^-4))
         maxx_tmp=maxx+flag_maxx*buchangx;
         kytest = 2*pi*(maxx_tmp-n)/(2*n);
         kxtest = 2*pi*(maxy-n)/(2*n);
         Irtest(:,:)=exp(1i*(kxtest*xx2+kytest*yy2));
         replcHtest(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1))).*Irtest(:,:))));
         replcHtest(abs(replcHtest)>0.9)=1;
         replcHtest(abs(replcHtest)~=1)=0;
         replch(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H))).*Irtest(:,:))));
         replch= replch.*replcHtest;
         replctest(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(spzhongyiwei))).*Irtest(:,:))));
         youhuatest=replctest.*replcHtest.*H;
         hetest=(conj(youhuatest)).*(spzhongxin.*replch);
         hetest=abs(sum(hetest(:)));
         cishuhtest=H1.*replcHtest;
         cishuhtest=sum(cishuhtest(:));
         hetest=hetest./cishuhtest;
         he_tmp=hetest;
         if(he_tmp<=he)
             buchangx=0.5*buchangx;
         elseif(he_tmp>he)                 
             he=he_tmp;
             maxx=maxx_tmp;
             break;
         end
     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%maxy%%%%%%%%%%%%%%%%%%%%
             maxy_tmp1=maxy-10^-5;
             maxy_tmp2=maxy+10^-5;
             kx_tmp1 = 2*pi*(maxy_tmp1-n)/(2*n);
             kx_tmp2 = 2*pi*(maxy_tmp2-n)/(2*n);
             for ii=1:1:2
                 switch ii
                     case 1
                         kxtest=kx_tmp1;
                         kytest=2*pi*(maxx-n)/(2*n);
                     case 2
                         kxtest=kx_tmp2;
                         kytest=2*pi*(maxx-n)/(2*n);
                 end
                 Irtest(:,:)=exp(1i*(kxtest*xx2+kytest*yy2));
                 replcHtest(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1))).*Irtest(:,:))));
                 replcHtest(abs(replcHtest)>0.9)=1;
                 replcHtest(abs(replcHtest)~=1)=0;
                 replch(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H))).*Irtest(:,:))));
                 replch= replch.*replcHtest;
                 replctest(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(spzhongyiwei))).*Irtest(:,:))));
                 youhuatest=replctest.*replcHtest.*H;
                 hetest=(conj(youhuatest)).*(spzhongxin.*replch);
                 hetest=abs(sum(hetest(:)));
                 cishuhtest=H1.*replcHtest;
                 cishuhtest=sum(cishuhtest(:));
                 hetest=hetest./cishuhtest;
                 test(ii)=hetest;
             end
             if((test(1)>test(2)))
                 flag_maxy=-1;
             elseif((test(1)<test(2)))
                 flag_maxy=+1;
             else
                 %error('Can not estimate the pattern wave vector');
                 flag_maxy=-1*flag_maxy;
             end
    while(buchangy>(10^-4))
        maxy_tmp=maxy+flag_maxy*buchangy;
        kytest = 2*pi*(maxx-n)/(2*n);
        kxtest = 2*pi*(maxy_tmp-n)/(2*n);
        Irtest(:,:)=exp(1i*(kxtest*xx2+kytest*yy2));
        replcHtest(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H1))).*Irtest(:,:))));
        replcHtest(abs(replcHtest)>0.9)=1;
        replcHtest(abs(replcHtest)~=1)=0;
        replch(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(H))).*Irtest(:,:))));
        replch= replch.*replcHtest;
        replctest(:,:) = fftshift(fft2(ifftshift(fftshift(ifft2(ifftshift(spzhongyiwei))).*Irtest(:,:))));
        youhuatest=replctest.*replcHtest.*H;
        hetest=(conj(youhuatest)).*(spzhongxin.*replch);
        hetest=abs(sum(hetest(:)));
        cishuhtest=H1.*replcHtest;
        cishuhtest=sum(cishuhtest(:));
        hetest=hetest./cishuhtest;
        he_tmp=hetest;
        if(he_tmp<=he)
            buchangy=0.5*buchangy;
        elseif(he_tmp>he)
            he=he_tmp;
            maxy=maxy_tmp;
            break;
        end
    end
end
 zuobiaox(spi+1+floor((spi-1)/2),:)=maxx;
 zuobiaoy(spi+1+floor((spi-1)/2),:)=maxy;
 zuobiaox(spi+2+floor((spi-1)/2),:)=2*n-maxx;
 zuobiaoy(spi+2+floor((spi-1)/2),:)=2*n-maxy;
 replc6(:,:,spi)=lihe;
 waitbar(spi/((nangles)*(nphases-1)), Progressbar, 'Parameter Estimate...');
end


 %save zuobiaox zuobiaox;
 %save zuobiaoy zuobiaoy;
 save([pathname  filename(1:end-4) '_zuobiaox.mat'],'zuobiaox');
 save([pathname  filename(1:end-4) '_zuobiaoy.mat'],'zuobiaoy');
%  imwritestack(abs(replc6),'lihe.tif');
 waitbar(1, Progressbar, 'Done');
 clearvars -except wenjianming pathname filename psf n weilac spjg bgflag bgname starframe beishu_an beishu_re filename_all filei  wavelengh pg fanwei regul deconv_otfname zstack testreadx testready savec notch_swith notch_para_1 notch_para_2  sigma mu nphases nangles Progressbar
 flag=1;
 SIM_512_fast_22_wl9;
  

