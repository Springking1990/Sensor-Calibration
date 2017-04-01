close all;
clc;
clear;
%���㷨������ϵ������Ҫ��һ��K,��ʾһ�������㵽��K�������㣻������Ⱥ��ռ������ı���Percent
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

RPtnum = 100;  % the number of random pointers ȫ����Ⱥ������ �û�ָ��
%************ generate random points ********����ģ��ȫ����Ⱥ��
minval = min(min(magn));
maxval = max(max(magn));

randPt = randi([minval maxval],RPtnum,3);
magn = [magn;randPt];
%*********** generate random points ********����ģ��ȫ����Ⱥ��

row = size(magn,1);
x = magn(:,1);
y = magn(:,2);
z = magn(:,3);

% the Kth nearest neigbour K>=2��K=1 represents distance between data and itself,in fact the distance=0��
% K=2 means the distance between data and its closet data
K = 10; %in fact, real k = K-1  �û�ָ��

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
% �ҳ�ÿ���������K������ھ��룬���������������
[Dk,Ind]=sort(MagSort(K,:),'descend');
%���þ�����ֵ �������ðٷֱ�
Dthreshold = 4.5e+4;

Percent = 0.03;%��Ϊ������ǰPercent*100%������Ⱥ�㡣  �û�ָ��
num=ceil(length(Ind)*Percent); % ceil()����ȡ��
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
magn(DeleOutl,:) = [];  %ɾ��ԭ�����б���������Ⱥ��
figure
plot3(magn(:,1),magn(:,2),magn(:,3),'.b');
hold on;
axis equal;
grid on;



