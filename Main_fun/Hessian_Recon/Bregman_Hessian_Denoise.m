%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 7. Code Structure: Scripts vs. functions 
% Make a script into a function for the MATLAB runtime to explore a reduced search space.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Bregman_Hessian_Denoise(mu, sigma, pathname, filename)
warning off
mkdir([pathname '\SIM-Hessian']);
warning on
%mu=str2double(get(handles.edit5,'String'));
%sigma=str2double(get(handles.edit6,'String'));
%%
clearvars -except filename pathname sigma mu
disp('Hessian Denoise,please wait...');
Progressbar = waitbar(0, 'Hessian reconstruction');
filename_notif=filename(1:end-4);
y = imreadstack([pathname filename]); %observed data 
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

labmda_div = gpuArray(single((1/lamda)));
s_lambda = gpuArray(single(siranu/lamda));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 6. Floating point number precision (single precision)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
y = gpuArray(single(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 5. Duplicated operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Load if a file exists
%else calcuate and store
diff_path = strcat(pathname, 'SIM-Hessian\FFT_diff\');
filename = strcat(diff_path, num2str(sx), '-', num2str(sy), '-', num2str(sz),'.mat');

if exist(filename, 'file') == 2
    %tic
    load(filename, 'divide');
    divide = gpuArray(divide);
else
    ztiduzz(:,:,1)=1;
    ztiduzz(:,:,2)=-2;
    ztiduzz(:,:,3)=1;

    ztiduxz(:,:,1)=[1,-1];
    ztiduxz(:,:,2)=[-1,1];

    ztiduyz(:,:,1)=[1;-1];
    ztiduyz(:,:,2)=[-1;1]; 

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
    divide = gpuArray(single((siranu/lamda) + Frefft));
    
    if exist(diff_path, 'dir') == 0
        mkdir(diff_path);
    end

    save(filename, 'divide');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 6. Floating point number precision (single precision)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% iteration
b1 = gpuArray(zeros(sizex,'single'));
b2 = b1;
b3 = b1;
b4 = b1;
b5 = b1;
b6 = b1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 5. Duplicated operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% The matrix is repeatedly computed in the loop in the original code. 
% A variable "frac_ori" stores the original value.
frac_ori = gpuArray((s_lambda)*(y)); 
frac = frac_ori;

% The values of diff_x11 and diff_x12 are used multiple times in the loop. 
% Variables "diff_x11" and "diff_x12" store the original values
diff_x11 = gpuArray(zeros(sizex,'single'));
diff_x12 = diff_x11;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 3.Sub. In-place memory operation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To avoid creating a matrix at runtime, all the functions use 
% pre-allocated matrixes, temp1 and temp2
temp1 = gpuArray(zeros(sizex, 'single'));
temp2 = temp1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Identifying algorithm performance bottlenecks.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the loop that incurs performance bottlenecks due to
% compute-intensive and iterative.
% Directions (xx, xy, xz, yy, yz, and zz) can be executed in CPU in parallel
% in separate functions as discussed in Section 8.
% In GPU usage (Section 9), all GPU cores are utilized without using
% separate functions as each task is executed asynchronously.
for ii = 1:iter_Bregman
        frac = fftn(frac);
       
        if ii>1
            x = real(ifftn(frac./divide));
            
            if ii == iter_Bregman
                break;
            end
        else
            x = real(ifftn(frac./(s_lambda)));
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Section 3.Sub. In-place memory operation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % In the loop, pre-allocated matrixes diff_x1(*) and temp(*)
    % are continuously used, avoding memory allocation during execution.
    diff_x11 = diff_x11.*0;
    diff_x12 = diff_x12.*0;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 3. Memory access and built-in MATLAB functions 
    % 4. Inline code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Instead of calling back_diff and forward_diff
	% Frequent function calls will degrade performance due to memory access (stack)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % XY direction
    diff_x11(1:sx-1,:,:) = diff(x, 1, 1); 
	temp2 = temp2.*0;
    temp2(2:sx,:,:) = diff(diff_x11, 1, 1);
    
    tb = temp2+b1;
    signd = abs(tb)-(labmda_div);
    signd = max(signd, 0);
    signd=signd.*sign(tb);
    b1 = b1+(temp2-signd);
    
    temp1 = temp1.*0;
    temp1(1:sx-1,:,:) = diff((signd-b1), 1, 1);
	temp2 = temp2.*0;
    temp2(2:sx,:,:) = diff(temp1, 1, 1);
    frac = frac_ori + temp2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% XY direction
	diff_x12(:,1:sy-1,:) = diff(x, 1, 2);
    temp2 = temp2.*0;
    temp2(:,2:sy,:) = diff(diff_x12, 1, 2);
    
    tb = temp2+b2;
    signd = abs(tb)-(labmda_div);
    signd = max(signd, 0);
    signd=signd.*sign(tb);
    b2 = b2+(temp2-signd);
        
    temp1 = temp1.*0;
    temp2 = temp2.*0;
    
    temp1(:,1:sy-1,:) = diff((signd-b2), 1, 2);
    temp2(:,2:sy,:) = diff(temp1, 1, 2);
    frac = frac + temp2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%XZ direction
    temp1 = temp1.*0;
    temp2 = temp2.*0;
    
    temp1(:,:,1:sz-1) = diff((x), 1, 3);
    temp2(:,:,2:sz) = diff(temp1, 1, 3);
        
    tb = temp2+b3;
    signd = abs(tb)-(labmda_div);
    signd = max(signd, 0);
    signd=signd.*sign(tb);
    b3 = b3+(temp2-signd);
    
    temp1 = temp1.*0;
    temp2 = temp2.*0;
    temp1(:,:,1:sz-1) = diff((signd-b3), 1, 3);
    temp2(:,:,2:sz) = diff(temp1, 1, 3);
    frac = frac + (zbei^2)*temp2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%%%%%%%%%
	%YY direction
    temp2 = temp2.*0;
    temp2(:,1:sy-1,:) = diff(diff_x11, 1, 2);
    
    tb = temp2 + b4;
    signd = abs(tb)-(labmda_div);
    signd = max(signd, 0);
    signd=signd.*sign(tb);
    b4 = b4+(temp2-signd);
    
    temp1 = temp1.*0;
    temp2 = temp2.*0;
    temp1(:,2:sy,:) = diff((signd-b4), 1, 2);
    temp2(2:sx,:,:) = diff(temp1, 1, 1);
    frac = frac + 2 * temp2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%YZ direction
    temp2 = temp2.*0;
    temp2(:,:,1:sz-1) = diff(diff_x11, 1, 3);
    
    tb = temp2+b5;
    signd = abs(tb)-(labmda_div);
    signd = max(signd, 0);
    signd=signd.*sign(tb);
    b5 = b5+(temp2-signd);
 
    temp1 = temp1.*0;
    temp2 = temp2.*0;
    temp1(:,:,2:sz) = diff((signd-b5), 1, 3);
    temp2(2:sx,:,:) = diff(temp1, 1, 1);
    frac = frac + 2 * (zbei)* temp2;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%ZZ direction
    temp2 = temp2.*0;
    temp2(:,:,1:sz-1) = diff(diff_x12, 1, 3);
    
    tb = temp2+b6;
    signd = abs(tb)-(labmda_div);
    signd = max(signd, 0);
    signd=signd.*sign(tb);
    b6 = b6+(temp2-signd);
    
    temp1 = temp1.*0;
    temp2 = temp2.*0;
    temp1(:,:,2:sz) = diff((signd-b6), 1, 3);
    temp2(:,2:sy,:) = diff(temp1, 1, 2);
    frac = frac + 2 * (zbei) * temp2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

    waitbar(ii/iter_Bregman , Progressbar, 'Hessian reconstruction');
end

toc
x(x<0) = 0;
x=x(:,:,1:y_flag);

% gather function will move data in GPU memory into CPU memory. 
x=gather(x);
imwritestack(x.*ymax, [pathname '\SIM-Hessian\Hessian-Denoise-' filename_notif '.tif']);
disp('Hessian reconstruction Successfully');
close(Progressbar);
helpdlg('All Done');