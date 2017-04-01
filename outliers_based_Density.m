close all;
clc;
clear;
%*****************�����ܶȵ���Ⱥ�����㷨*******************
%�����û�ָ���Ĳ���ΪK��
%�ֲ���Ⱥ������LoFk����ֵ�Ͻ�Big_threshold����ֵ�½�Small_threshold
%�ֲ��ɴ��ܶ�Lrdk ����Percent

%*****************�����ܶȵ���Ⱥ�����㷨*******************

%magn = load('magn_rawdata1.txt');
%magn = load('magn_rawdata2.txt');
%magn = load('magn_rawdata3.txt');
magn = load('magn_rawdata4.txt');
%magn = load('magn_rawdata5.txt');
%magn = load('magn_rawdata6.txt');
%magn = load('magn_rawdata7.txt');
%magn = load('magn_rawdata8.txt');
%magn = load('magn_rawdata9.txt');
%magn = load('magn_result_8.txt');
%magn = load('magn_rawdata10.txt');
%magn = load('magn_result_11.txt');
%magn = load('magn_result_C.txt');
%magn = load('nonmagn_result_1.txt');

RPtnum = 100;  % the number of random pointers ȫ����Ⱥ������  �û�ָ��
%************ generate random points ******����ģ��ȫ����Ⱥ��
minval = min(min(magn));
maxval = max(max(magn));

randPt = randi([minval maxval],RPtnum,3);
magn = [magn;randPt];
%*********** generate random points *******����ģ��ȫ����Ⱥ��

row = size(magn,1);
x = magn(:,1);
y = magn(:,2);
z = magn(:,3);

% the Kth nearest neigbour K>=2��K=1 represents distance between data and itself,in fact the distance=0��
% K=2 means the distance between data and its closet data
K = 20;  %in fact, real k = K-1   �û�ָ��

%Distmat=zeros(row,row); % the matrix is applied to computer
mag_x  = repmat(magn(:,1),1,row);
mag_y  = repmat(magn(:,2),1,row);
mag_z  = repmat(magn(:,3),1,row);
mag_tx  = mag_x';
mag_ty  = mag_y';
mag_tz  = mag_z';
%distmat(i,j) represents the eulerdistance between point i and point j;
Distmat = sqrt((mag_x-mag_tx).^2+(mag_y-mag_ty).^2+(mag_z-mag_tz).^2); 

% magSort represents sorted elements in each column by ascending order
%index is the same size as distmat and describes the arrangement of the elements of distmat into magSort along the sorted dimension
[MagSort,Index]=sort(Distmat,'ascend');

% Nk represents sets of kth nearest neigbor point,Nk_num stores distances
% Nk_ind stored index of each match
Nk_num = MagSort(2:K,:);
Nk_ind = Index(2:K,:);
dic_y = 1;
dic_x = K-1;
% �ҳ�����o��k����ڵļ��ϣ���Ϊ�����ж������o�ľ�����ȣ�������Ҫ�жϣ�
while(dic_y<=size(magn,1))
    if(Nk_num(dic_x,dic_y)==MagSort(dic_x+2,dic_y))
        Nk_num(dic_x+1,dic_y) = MagSort(dic_x+2,dic_y);
        Nk_ind(dic_x+1,dic_y) = Index(dic_x+2,dic_y);
        dic_x = dic_x+1;
    else
        dic_y = dic_y+1;
        dic_x = K-1;
    end
end

% ����o�ɴ����reachdisk(o'-o)  Nk_cntΪ||Nk(o)||
distK=MagSort(K,:);% ÿ���������k���������
Nk_column = 1;
Nk_row = 1;
reachdist=zeros(length(Nk_ind(:,1)),length(Nk_ind(1,:)));%Ԥ�����ڴ棬����������ٶ�
Nk_cnt = zeros(1,row);
while(Nk_column <= size(magn,1))
    if((Nk_row <= length(Nk_ind(:,1))) && (Nk_ind(Nk_row,Nk_column) ~= 0))
        reachdist(Nk_row,Nk_column) = max(distK(Nk_ind(Nk_row,Nk_column)),Nk_num(Nk_row,Nk_column));
        %reachdist(Nk_row,Nk_column) = Nk_num(Nk_row,Nk_column);
        Nk_row = Nk_row+1;
    else
        Nk_cnt(1,Nk_column) = length(find(Nk_ind(:,Nk_column) ~= 0));%ͳ����������k�ڽ��򼯺���Ԫ�ظ���
        Nk_column = Nk_column+1;
        Nk_row = 1;
    end       
end
% �������o�ľֲ��ɴ��ܶ�lrdk(o)
Init(1,1:length(reachdist(:,1)))=1;
Sigmard = Init*reachdist;
lrdk = Nk_cnt/(diag(Sigmard));

% �������o�ľֲ���Ⱥ������LOFk(o)��
Imat=zeros(row,row);%Ԥ�����ڴ棬����������ٶ�
Imat_row=1;
count=1;
while(Imat_row<=row)
    if(count<=Nk_cnt(1,Imat_row))
       Imat(Imat_row,Nk_ind(count,Imat_row))=1;
       count=count+1; 
    else
        count = 1;
        Imat_row = Imat_row+1;
    end      
end
sigma_lrdk = Imat*lrdk';
LOFk=(diag(Nk_cnt))^2 \(diag(Sigmard)*(sigma_lrdk));

[Den_num,Den_ind]=sort(LOFk,'descend');

%�Ծֲ���Ⱥ������LOFk��Ϊ����ָ��  �ʺϼ�⡾�ܼ������д��ڵ�ϡ����Ⱥ�㡿
%���þֲ���Ⱥ��������ֵ,����������Ⱥ����ռ�ٷֱ�
Big_threshold = 1.5;    %������ֵ�Ͻ�   �û�ָ��
smal_threshold = 0.75;    %������ֵ�½�    �û�ָ��
Den_row=1;
figure
plot3(magn(:,1),magn(:,2),magn(:,3),'.b');
hold on;
axis equal;
grid on;
Outliers = 0;
for paint_cnt=1:row
    if(Den_num(paint_cnt)>=Big_threshold)
        Outliers = [Outliers Den_ind(paint_cnt)];  %LOFk���Ӽ�������Ⱥ��
        plot3(x(Den_ind(paint_cnt)),y(Den_ind(paint_cnt)),z(Den_ind(paint_cnt)),'rp');
    else
        if(Den_num(paint_cnt)<=smal_threshold)
            plot3(x(Den_ind(paint_cnt)),y(Den_ind(paint_cnt)),z(Den_ind(paint_cnt)),'go');
        end
    end
    hold on;
end
%{
while(Den_num(Den_row)>=Big_threshold)
    plot3(x(Den_ind(Den_row)),y(Den_ind(Den_row)),z(Den_ind(Den_row)),'rp');
    Den_row=Den_row+1;
    hold on;
end
%}

%�Զ���o�ľֲ��ɴ��ܶ�lrdk��Ϊ����ָ��   �ʺϼ�⡾�ܶȸ߻�͵���Ⱥ�㡿
[lrdk_num,lrdk_ind] = sort(lrdk,'descend');

Percent = 0.05;%��Ϊ������ǰPercent*100%������Ⱥ�㡣  �û�ָ��
lrdk_cnt=ceil(length(lrdk_ind)*Percent); % ceil()����ȡ��

figure
plot3(magn(:,1),magn(:,2),magn(:,3),'.b');
hold on;
axis equal;
grid on;
for i=1:lrdk_cnt
    plot3(x(lrdk_ind(i)),y(lrdk_ind(i)),z(lrdk_ind(i)),'ro');
    hold on;
end

DeleOutl = nonzeros(unique([Outliers lrdk_ind(1:lrdk_cnt)]));
magn(DeleOutl,:) = [];  %ɾ��ԭ�����б���������Ⱥ��
figure
plot3(magn(:,1),magn(:,2),magn(:,3),'.b');
hold on;
axis equal;
grid on;











