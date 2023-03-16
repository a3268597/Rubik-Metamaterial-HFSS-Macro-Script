%Direction. For metasurface with vortex.
%%Feb. first , 2023. by Xiaokun Yang
clear, clc
units = 'mm';
Units = 'mm';
Type = 'csv';

%%%%%%%Initialization
fc=18;
M=20;
N=M;
arraynum=M*N;
%design Variables
    f=18;  %GHz
    pp=10; %period
    b2=pp;
    da=12;
%%此段为计算不同模式阵列相位分布(squar_loop)。a矩阵为输出相位矩阵

a=ones(N,M);%a矩阵为输出相位矩阵
focal_diameter_ratio = 0.8;%焦径比的系数  从0.2-1；
frequency = 18;%中心频率
P=10e-3;%周期
compensation_horn_distance = 76.2e-3;
L=-1;%算高增益时L为0
arraynum=M*N;
phi=0;
theta1=0;
%变量的值
%%%%%%%%%%对生成的矩阵每个位置所需相位进行计算
ArrayP =  P*(M/2)-P/2;
x1 = -ArrayP : P : ArrayP;
y1 = -ArrayP : P : ArrayP;
z2 = focal_diameter_ratio * sqrt((P*M)^2+(P*M)^2);
wavelength = 0.3/frequency;
k0 = (2*pi/wavelength);

for i =1:M
    for j=1:N        
        %A side
        calculate_phaseA(i,j)=k0*(sqrt(x1(i)^2+y1(j)^2+z2^2)-z2)+0*atan2(y1(j),x1(i));
        calculate_phaseA(i,j)=calculate_phaseA(i,j)*180/pi;
        calculate_phaseA(i,j)=mod(calculate_phaseA(i,j),360);
        %B side
        calculate_phaseB(i,j)=k0*(sqrt(x1(i)^2+y1(j)^2+z2^2)-z2)+1*atan2(y1(j),x1(i));
        calculate_phaseB(i,j)=calculate_phaseB(i,j)*180/pi;
        calculate_phaseB(i,j)=mod(calculate_phaseB(i,j),360);
        %C side
        calculate_phaseC(i,j)=k0*(sqrt(x1(i)^2+y1(j)^2+z2^2)-z2)+2*atan2(y1(j),x1(i));
        calculate_phaseC(i,j)=calculate_phaseC(i,j)*180/pi;
        calculate_phaseC(i,j)=mod(calculate_phaseC(i,j),360);
        %D side
        calculate_phaseD(i,j)=k0*(sqrt(x1(i)^2+y1(j)^2+z2^2)-z2)+3*atan2(y1(j),x1(i));
        calculate_phaseD(i,j)=calculate_phaseD(i,j)*180/pi;
        calculate_phaseD(i,j)=mod(calculate_phaseD(i,j),360);
        %E side
        calculate_phaseE(i,j)=k0*(sqrt(x1(i)^2+y1(j)^2+z2^2)-z2)+(-1)*atan2(y1(j),x1(i));
        calculate_phaseE(i,j)=calculate_phaseE(i,j)*180/pi;
        calculate_phaseE(i,j)=mod(calculate_phaseE(i,j),360);
        %F side
        calculate_phaseF(i,j)=k0*(sqrt(x1(i)^2+y1(j)^2+z2^2)-z2)+(-2)*atan2(y1(j),x1(i));
        calculate_phaseF(i,j)=calculate_phaseF(i,j)*180/pi;
        calculate_phaseF(i,j)=mod(calculate_phaseF(i,j),360);

    end
end

%Mapping
p1=0.02313;
p2=3.006;
for i =1:N
    for j=1:M
        calculate_phaseA(i,j)=p1*calculate_phaseA(i,j)+p2;
        calculate_phaseB(i,j)=p1*calculate_phaseB(i,j)+p2;   
        calculate_phaseC(i,j)=p1*calculate_phaseC(i,j)+p2;   
        calculate_phaseD(i,j)=p1*calculate_phaseD(i,j)+p2;   
        calculate_phaseE(i,j)=p1*calculate_phaseE(i,j)+p2;   
        calculate_phaseF(i,j)=p1*calculate_phaseF(i,j)+p2;   
    end
end

%%%%%%%%%脚本开始准备
false = 0;
true = 1;
fileName = ['squareloopN',num2str(N),'_',num2str(fc),'G_',num2str(arraynum),'N_l=',num2str(6)];
tmpScriptFile = ['.\',fileName,'.vbs'];
%开始写VBS脚本
    fid = fopen(tmpScriptFile, 'wt');   % 'wt'表示以文本模式打开文件，可写，覆盖原有内容
   % 创建一个新的工程并插入一个新的设计
    hfssNewProject(fid);
    hfssInsertDesign(fid, fileName);
%     hfssPre(fid);
   
    hfssaddVar(fid, 'f0', f, []);
    hfssaddVar(fid, 'pp', pp, units);
    hfssaddVar(fid, 'b2', b2, units);
    hfssaddVar(fid, 'M', M, []);  
    hfssaddVar(fid, 'N', N, []); 
    %%%%%%%%Array of metasrface(square loop)
    squareloopNameA = cell(arraynum, 1);
    squareloopNameB = cell(arraynum, 1);
    squareloopNameC = cell(arraynum, 1);
    squareloopNameD = cell(arraynum, 1);
    squareloopNameE = cell(arraynum, 1);
    squareloopNameF = cell(arraynum, 1);
    for n=1:N
        for m=1:M
            i=N*(m-1)+n;
            %%%%add local CS
            CSName = ['CS', num2str(i)];
            Origin = {['(',num2str(m),' -(M+2)/2)*pp'],['((N)/2-',num2str(n),')*pp'],'0'};
            XAxisVec = [1 0 0 ];
            YAxisVec = {0, 1, 0};
            hfssCreateRelativeCS(fid, CSName, Origin,XAxisVec, YAxisVec, units);
            hfssSetWCS(fid, CSName);
   %%%%%%%Array of metasrface(square loop)
           %Side A
            squareloopNameA{i} = ['squareloopA_',num2str(i)];
            Start_squareloopA = {'0','0','M/2*pp'};
            Start_sizeA = {'b2','b2',calculate_phaseA(n,m)};
            hfssBox(fid, squareloopNameA{i}, Start_squareloopA, Start_sizeA, units); 
            hfssSetWCS(fid, 'Global');
       %Side B
            squareloopNameB{i} = ['squareloopB_',num2str(i)];
            Start_squareloopB = {'M/2*pp',(m-(M/2+1))*pp,(n-(M/2+1))*pp};
            Start_sizeB = {calculate_phaseB(n,m),'b2','b2'};
            hfssBox(fid, squareloopNameB{i}, Start_squareloopB, Start_sizeB, units); 
       %Side C
            squareloopNameC{i} = ['squareloopC_',num2str(i)];
            Start_squareloopC = {(m-(M/2+1))*pp,(n-(M/2+1))*pp,'-M/2*pp'};
            Start_sizeC = {'b2','b2',-calculate_phaseC(n,m)};
            hfssBox(fid, squareloopNameC{i}, Start_squareloopC, Start_sizeC, units); 
       %Side D
            squareloopNameD{i} = ['squareloopD_',num2str(i)];
            Start_squareloopD = {'-M/2*pp',(m-(M/2+1))*pp,(n-(M/2+1))*pp};
            Start_sizeD = {-calculate_phaseD(n,m),'b2','b2'};
            hfssBox(fid, squareloopNameD{i}, Start_squareloopD, Start_sizeD, units); 
       %Side E
           squareloopNameE{i} = ['squareloopE_',num2str(i)];
            Start_squareloopE = {(m-(M/2+1))*pp,'M/2*pp',(n-(M/2+1))*pp};
            Start_sizeE = {'b2',calculate_phaseE(n,m),'b2'};
            hfssBox(fid, squareloopNameE{i}, Start_squareloopE, Start_sizeE, units); 
       %Side F
            squareloopNameF{i} = ['squareloopF_',num2str(i)];
            Start_squareloopF = {(m-(M/2+1))*pp,'-M/2*pp',(n-(M/2+1))*pp};
            Start_sizeF = {'b2',-calculate_phaseE(n,m),'b2'};
            hfssBox(fid, squareloopNameF{i}, Start_squareloopF, Start_sizeF, units); 
        end
    end
fclose(fid);

