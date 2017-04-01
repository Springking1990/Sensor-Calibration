close all;
clc;
clear;
%该算法有两个系数很重要，一是K,表示一个样本点到第K个样本点；二是离群点占样本点的比例Percent
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

RPtnum = 100;  % the number of random pointers 全局离群点数量 用户指定
%************ generate random points ********用于模拟全局离群点
minval = min(min(magn));
maxval = max(max(magn));

randPt = randi([minval maxval],RPtnum,3);
magn = [magn;randPt];
%*********** generate random points ********用于模拟全局离群点

row = size(magn,1);
x = magn(:,1);
y = magn(:,2);
z = magn(:,3);

% the Kth nearest neigbour K>=2（K=1 represents distance between data and itself,in fact the distance=0）
% K=2 means the distance between data and its closet data
K = 10; %in fact, real k = K-1  用户指定

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
%MagComp=MagSort(1:K,:); 
%Indcomp=Index(1:K,:);
% 找出每个样本点第K个最近邻距离，并按降序进行排序
[Dk,Ind]=sort(MagSort(K,:),'descend');
%设置距离阈值 或者设置百分比
Dthreshold = 4.5e+4;

Percent = 0.03;%认为最大距离前Percent*100%都是离群点。  用户指定
num=ceil(length(Ind)*Percent); % ceil()向上取整
figure
plot3(magn(:,1),magn(:,2),magn(:,3),'.b');
hold on;
axis equal;
grid on;
for i=1:num
    plot3(x(Ind(i)),y(Ind(i)),z(Ind(i)),'rp');
    hold on;
end

DeleOutl = Ind(1:num);
magn(DeleOutl,:) = [];  %删除原样本中被检测出的离群点
figure
plot3(magn(:,1),magn(:,2),magn(:,3),'.b');
hold on;
axis equal;
grid on;



