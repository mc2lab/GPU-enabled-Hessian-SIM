%%read images
disp('Start reconstruction,please wait...');
pathname=handles.pathname;
filename=handles.filename;
warning off
mkdir([pathname '\Pseudo-TIRF']);
mkdir([pathname '\SIM-Wiener']);
mkdir([pathname '\SIM-Hessian']);
mkdir([pathname '\SIM-TV']);
mkdir([pathname '\Running-Average']);
warning on
avger(pathname,filename)
%%
% parameter
weilac=str2double(get(handles.edit1,'String'));
wavelengh=str2double(get(handles.edit2,'String'));
beishu_an=str2double(get(handles.edit3,'String'));
mu=str2double(get(handles.edit5,'String'));
sigma=str2double(get(handles.edit6,'String'));
Pixel_Size=str2double(get(handles.edit7,'String'));
Excition_NA=str2double(get(handles.edit8,'String'));
spjg=[str2double(get(handles.edit9,'String')),str2double(get(handles.edit10,'String')),str2double(get(handles.edit11,'String'))];
otf_flag=get(handles.radiobutton4,'Value')-get(handles.radiobutton5,'Value');
if get(handles.radiobutton6,'Value') && ~get(handles.radiobutton7,'Value')
    bgname='background.tif';
elseif ~get(handles.radiobutton6,'Value') && get(handles.radiobutton7,'Value')
    bgname=[handles.bg_pathname '\' handles.bg_filename];
end
if  get(handles.radiobutton1,'Value')==1 && get(handles.radiobutton2,'Value')==0 && get(handles.radiobutton3,'Value')==0
    nphases = 3; 
    nangles = 3; 
elseif  get(handles.radiobutton1,'Value')==0 && get(handles.radiobutton2,'Value')==0 && get(handles.radiobutton3,'Value')==1
    nphases = 3; 
    nangles = 2; 
end
beishu_re=1;
bgflag=1;
notch_swith=0;
notch_para_1=0.05;      
notch_para_2=1.2;       
savec=0;                %0 means not save the parameter c and angle£¬1 means save 
starframe=1;
%%
pg=512*2*Excition_NA*Pixel_Size/wavelengh;
info = imfinfo((fullfile(pathname, filename)));
zstack = numel(info);
testreadx = info(1).Height;
testready = info(1).Width;
clear  info;
testn = max([testreadx,testready,256]);

if wavelengh==488
    if otf_flag==1
        deconv_otfname='488OTF_512.tif';
    elseif otf_flag==-1
        deconv_otfname=[ handles.otf_pathname '\' handles.otf_filename ];
    else
        deconv_otfname='488OTF_512.tif';
    end
    psf=imreadstack(deconv_otfname);
    fanwei=50;
    regul=2*pi;
elseif wavelengh==561
    if otf_flag==1
        deconv_otfname='561OTF_512.tif';
    elseif otf_flag==-1
        deconv_otfname=[ handles.otf_pathname '\' handles.otf_filename ];
    else
        deconv_otfname='561OTF_512.tif';
    end
    psf=imreadstack(deconv_otfname);
    fanwei=50;
    regul=2*pi;
elseif wavelengh==647
    if otf_flag==1
        deconv_otfname='647OTF_512.tif';
    elseif otf_flag==-1
        deconv_otfname=[ handles.otf_pathname '\' handles.otf_filename ];
    else
        deconv_otfname='647OTF_512.tif';
    end
    psf=imreadstack(deconv_otfname);
    fanwei=50;
    regul=2*pi;
else
    if otf_flag==-1
        deconv_otfname=[ handles.otf_pathname '\' handles.otf_filename ];
    else
        warndlg('The default OTF only support 488/561/647 wavelength. Please choose your special OTF if other wavelength.','warn','modal'); 
        return;
    end
    psf=imreadstack(deconv_otfname);
    fanwei=50;
    regul=2*pi;
end
Progressbar = waitbar(0, 'Parameter Estimation...');
mysim3_512_fast_21