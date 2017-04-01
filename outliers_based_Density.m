close all;
clc;
clear;
%*****************基于密度的离群点检测算法*******************
%需用用户指定的参数为K，
%局部离群点因子LoFk的阈值上界Big_threshold、阈值下界Small_threshold
%局部可达密度Lrdk 参数Percent

%*****************基于密度的离群点检测算法*******************

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

RPtnum = 100;  % the number of random pointers 全局离群点数量  用户指定
%************ generate random points ******用于模拟全局离群点
minval = min(min(magn));
maxval = max(max(magn));

randPt = randi([minval maxval],RPtnum,3);
magn = [magn;randPt];
%*********** generate random points *******用于模拟全局离群点

row = size(magn,1);
x = magn(:,1);
y = magn(:,2);
z = magn(:,3);

% the Kth nearest neigbour K>=2（K=1 represents distance between data and itself,in fact the distance=0）
% K=2 means the distance between data and its closet data
K = 20;  %in fact, real k = K-1   用户指定

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
% 找出对象o的k最近邻的集合（因为可能有多个对象到o的距离相等，所以需要判断）
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

% 计算o可达距离reachdisk(o'-o)  Nk_cnt为||Nk(o)||
distK=MagSort(K,:);% 每个样本点的k最邻域距离
Nk_column = 1;
Nk_row = 1;
reachdist=zeros(length(Nk_ind(:,1)),length(Nk_ind(1,:)));%预分配内存，以提高运算速度
Nk_cnt = zeros(1,row);
while(Nk_column <= size(magn,1))
    if((Nk_row <= length(Nk_ind(:,1))) && (Nk_ind(Nk_row,Nk_column) ~= 0))
        reachdist(Nk_row,Nk_column) = max(distK(Nk_ind(Nk_row,Nk_column)),Nk_num(Nk_row,Nk_column));
        %reachdist(Nk_row,Nk_column) = Nk_num(Nk_row,Nk_column);
        Nk_row = Nk_row+1;
    else
        Nk_cnt(1,Nk_column) = length(find(Nk_ind(:,Nk_column) ~= 0));%统计样本对象k邻近域集合中元素个数
        Nk_column = Nk_column+1;
        Nk_row = 1;
    end       
end
% 计算对象o的局部可达密度lrdk(o)
Init(1,1:length(reachdist(:,1)))=1;
Sigmard = Init*reachdist;
lrdk = Nk_cnt/(diag(Sigmard));

% 计算对象o的局部离群点因子LOFk(o)。
Imat=zeros(row,row);%预分配内存，以提高运行速度
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

%以局部离群点因子LOFk作为评价指标  适合检测【密集数据中存在的稀疏离群点】
%设置局部离群点因子阈值,或者设置离群点所占百分比
Big_threshold = 1.5;    %设置阈值上界   用户指定
smal_threshold = 0.75;    %设置阈值下界    用户指定
Den_row=1;
figure
plot3(magn(:,1),magn(:,2),magn(:,3),'.b');
hold on;
axis equal;
grid on;
Outliers = 0;
for paint_cnt=1:row
    if(Den_num(paint_cnt)>=Big_threshold)
        Outliers = [Outliers Den_ind(paint_cnt)];  %LOFk因子检测出的离群点
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

%以对象o的局部可达密度lrdk作为评价指标   适合检测【密度高或低的离群点】
[lrdk_num,lrdk_ind] = sort(lrdk,'descend');

Percent = 0.05;%认为最大距离前Percent*100%都是离群点。  用户指定
lrdk_cnt=ceil(length(lrdk_ind)*Percent); % ceil()向上取整

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
magn(DeleOutl,:) = [];  %删除原样本中被检测出的离群点
figure
plot3(magn(:,1),magn(:,2),magn(:,3),'.b');
hold on;
axis equal;
grid on;











