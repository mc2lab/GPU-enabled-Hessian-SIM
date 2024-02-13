warning off stats:regress:NoConst
R2_a_ang=replc6_ang(:,:,ii);
R2_b_ang=R2_a_ang(:);
R2_b_ang(R2_b_ang~=0)=1;
R2_c=R2_ang(:,:,ii);
R2_d=R2_c(:);
R2_b_ang(:,2)=R2_d;
R2_c=R2_zhongxin_ang(:,:,ii);
R2_d=R2_c(:);
R2_b_ang(:,3)=R2_d;
R2_b_ang((R2_b_ang(:,1)==0),:)=[];
R2_y=R2_b_ang(:,2);
R2_x=R2_b_ang(:,3);
[B,bint,r,rint,stats]=regress(R2_y,R2_x);
R2_angles(ii)=stats(1);
if R2_angles(ii) < 0.1
    close(Progressbar);
    warndlg('The R-square statistic of linear regress of parameters is lower than 0.1','warn','modal'); 
    return;
end
warning on stats:regress:NoConst
