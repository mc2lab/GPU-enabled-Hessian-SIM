%update from Bregman_Hessian_LowRam,add 2*xy
%%
pathname=handles.pathname;
filename=handles.filename;
warning off
mkdir([pathname '\SIM-Hessian']);
warning on
mu=str2double(get(handles.edit5,'String'));
sigma=str2double(get(handles.edit6,'String'));
%%
clearvars -except filename pathname sigma mu
disp('Hessian Denoise,please wait...');
Progressbar = waitbar(0, 'Hessian reconstruction');
filename_notif=filename(1:end-4);
y = (imreadstack([pathname filename])); %observed data 
%% initialization
iter_Bregman = 1e2;     %number of iteration
lamda = 0.5;              
gexiang=1;              
siranu=mu;
zbei=sigma;
tic
y_flag=size(y,3);
if y_flag<3
    zbei=0;
    y(:,:,end+1:end+(3-size(y,3)))=repmat(y(:,:,end),[1,1,3-size(y,3)]);
    msgbox('Number of data frame is smaller than 3, the t and z-axis of Hessian was turned off(sigma=0)'); 
end
ymax=max(y(:));
y=y./ymax;
[sx,sy,sz] = size(y);
sizex=[sx,sy,sz] ;
x = zeros(sizex);                  %start point

ztiduzz(:,:,1)=1;
ztiduzz(:,:,2)=-2;
ztiduzz(:,:,3)=1;

ztiduxz(:,:,1)=[1,-1];
ztiduxz(:,:,2)=[-1,1];

ztiduyz(:,:,1)=[1;-1];
ztiduyz(:,:,2)=[-1;1];

%FFT of difference operator
tmp_fft=fftn([1 -2 1],sizex).*conj(fftn([1 -2 1],sizex));
Frefft = tmp_fft;
tmp_fft=fftn([1 ;-2 ;1],sizex).*conj(fftn([1; -2 ;1],sizex));
Frefft=Frefft + tmp_fft;
tmp_fft=fftn(ztiduzz,sizex).*conj(fftn(ztiduzz,sizex));
Frefft=Frefft +(zbei^2)*tmp_fft;
tmp_fft=fftn([1 -1;-1 1],sizex).*conj(fftn([1 -1;-1 1],sizex));
Frefft=Frefft + 2 * tmp_fft;
tmp_fft=fftn(ztiduxz,sizex).*conj(fftn(ztiduxz,sizex));
Frefft=Frefft + 2 * (zbei)*tmp_fft;
tmp_fft= fftn(ztiduyz,sizex).*conj(fftn(ztiduyz,sizex));
Frefft=Frefft + 2 * (zbei)*tmp_fft;
clear  tmp_fft
divide = single((siranu/lamda) + Frefft);
clear  Frefft
%% iteration
b1 = zeros(sizex,'single');
b2 = zeros(sizex,'single');
b3 = zeros(sizex,'single');
b4 = zeros(sizex,'single');
b5 = zeros(sizex,'single');
b6 = zeros(sizex,'single');
x = zeros(sizex,'int32');
frac = (siranu/lamda)*(y); 
for ii = 1:iter_Bregman
%% renew x
        frac = fftn(frac);
        if ii>1
            x = real(ifftn(frac./divide));
        else
            x = real(ifftn(frac./(siranu/lamda)));
        end

%% calculate the dirivative of x
%% renew d
% 'gexiang == 1' means anisotropic;'otherwise' means isotropic      
    frac = (siranu/lamda)*(y); 
    u = back_diff(forward_diff(x,1,1),1,1);
    signd = abs(u+b1)-1/lamda;
    signd(signd<0)=0;
    signd=signd.*sign(u+b1);
    d=signd;
    b1 = b1+(u-d);
    frac = frac+back_diff(forward_diff(d-b1,1,1),1,1);

    u = back_diff(forward_diff(x,1,2),1,2);
    signd = abs(u+b2)-1/lamda;
    signd(signd<0)=0;
    signd=signd.*sign(u+b2);
    d=signd;
    b2 = b2+(u-d);
    frac = frac+back_diff(forward_diff(d-b2,1,2),1,2);

    u = back_diff(forward_diff(x,1,3),1,3);
    signd = abs(u+b3)-1/lamda;
    signd(signd<0)=0;
    signd=signd.*sign(u+b3);
    d=signd;
    b3 = b3+(u-d);
    frac = frac+(zbei^2)*back_diff(forward_diff(d-b3,1,3),1,3);

    u = forward_diff(forward_diff(x,1,1),1,2);
    signd = abs(u+b4)-1/lamda;
    signd(signd<0)=0;
    signd=signd.*sign(u+b4);
    d=signd;
    b4 = b4+(u-d);
    frac = frac+ 2 * back_diff(back_diff(d-b4,1,2),1,1);

    u = forward_diff(forward_diff(x,1,1),1,3);
    signd = abs(u+b5)-1/lamda;
    signd(signd<0)=0;
    signd=signd.*sign(u+b5);
    d=signd;
    b5 = b5+(u-d);
    frac = frac+ 2 * (zbei)*back_diff(back_diff(d-b5,1,3),1,1);

    u = forward_diff(forward_diff(x,1,2),1,3);
    signd = abs(u+b6)-1/lamda;
    signd(signd<0)=0;
    signd=signd.*sign(u+b6);
    d=signd;
    b6 = b6+(u-d);
    frac = frac+ 2 * (zbei)*back_diff(back_diff(d-b6,1,3),1,2);
%     ii
    waitbar(ii/iter_Bregman , Progressbar, 'Hessian reconstruction');
end
toc
x(x<0) = 0;
x=x(:,:,1:y_flag);
imwritestack(x.*ymax, [pathname '\SIM-Hessian\Hessian-Denoise-' filename_notif '.tif']);
disp('Hessian reconstruction Successfully');
close(Progressbar);
helpdlg('All Done');