% function [number_aveS,number_aveE,number_aveI,number_aveH,number_aveR]=I_f_r()  %函数输入为网络的邻接矩阵与连接扰动强度，输出稳态时各状态节点的比例
clc;clear;
%1.生成初始邻接矩阵
load('BAnetwork1.mat');  %导入WS网络及邻接矩阵   %此处改变导入的网络、邻接矩阵
A0=A;
% %2.参数设置
%S的三种转化状态
lambda1=[0,0.2,0.4];        %中立——S变为E的概率 lambda1
lambda2=0.5;       %支持——S变为I态的概率lambda2
lambda3=0.1;       %反对——S变为I态的概率 lambda3(f)
% %记忆遗忘机制（H与I之间的转化）
delta=0.2;       %遗忘——I变为H的概率  delta
xi=0.5;       %接触记忆——H被动变为I的概率 ksai
eta=0.5;       %自发记忆——H主动变为I的概率 eta
fa=0.4;  %phi
epsilon=0.4;    %epsilon
T=30;                  %T为谣言在网络中扩散的时间总长，时间步长为1
s_number=100;          %s为仿真次数，对最后结果，通过多次仿真求平均值
y=size(lambda1,2);%取lambda1的列数
number_basic_I=zeros(y,T);
%保存各次仿真结果
number_S=zeros(s_number,T);    
number_E=zeros(s_number,T);
number_I=zeros(s_number,T);
number_H=zeros(s_number,T);
number_R=zeros(s_number,T);
% number_N=zeros(s_number,T);
N=length(A0);        %初始化N
pp=0.2
%3.仿真过程
for f=1:y;
for z=1:s_number
    s=zeros(N,T);       %初始化s，s为一个二维数组，用于存放每次扩散后，网络中各个节点的状态，维数也需要动态变化
 %生成初始感染者---------------------------
 ch_q= randperm(N,N*pp);
 for i=1:N*pp;
   countadj=find( A0(ch_q(i),:)==1 );%找到邻居节点所在的列
   sumadj=length(countadj);%邻居节点的数量
   ch_q2=ceil(rand(1,1)*sumadj);%随机选择邻居节点的列
   sd=countadj(ch_q2);%返回免疫者所在的行
   s(sd,1)=4;
 end
 intial=find(s(:,1)==0);
 sumintal=length(intial);
 ch_q3=ceil(rand(1,1)*sumintal);%随机选择邻居节点的列
 sd2=intial(ch_q3);
 s(sd2,1)=2;
 countr= find(s(:,1)==4);
sumr=length(countr)./N;%显示被免疫的节点站总节点的比例
%4信息扩散过程
 %每经历一次扩散后，对网络是否达到稳态进行判断
 for t=1:T-1
%      numI=length(find(s(:,t)==2));       
%         if (numI==0)                         %如果t时刻感染节点的个数为0，就把t时刻以后所有节点状态置为与t时刻相同，即传播达到了稳态
%             for h=t+1:T
%                 s(:,h)=s(:,t);
%             end
%             break;
%         end
        for i=1:N
 %-------------健康节点以alpha被传染，lambda1-3决定其转移状态---------------------
             if  s(i,t)==0
                 k=0;                       %统计健康节点的邻居节点中是感染节点的数量
                 Adjac1=find( A0(i,:)==1 ); %找出节点i的邻居节点位置，就在这里邻接矩阵换成免疫后的新矩阵A0
                 num_Adjac1=length(Adjac1); %计算节点i的度
                 for g=1:num_Adjac1
                     if s(Adjac1(g),t)==2    %%%%
                        k=k+1;             %统计健康节点的邻居节点中是感染节点的数量k
                     end
                 end
                p1=rand(1,1);
                p2=rand(1,1);
                p3=rand(1,1);
                s1=1-(1-lambda1(f))^k;
                s2=1-(1-lambda2)^k;
                s3=1-(1-lambda3)^k;
                if p1<=s1
                    s(i,t+1)=1;  %S变为E
                elseif p2<=s2
                    s(i,t+1)=2;  %S变为I态
                elseif p3<=s3;
                    s(i,t+1)=4;  %S变为R态
                end
             end
%----------------犹豫机制-E以phi直接变为I态(无需状态间连接)--------------------
            if  s(i,t)==1
                p=rand(1,1);
                if p<=fa
                    s(i,t+1)=2;
                else s(i,t+1)=1;
                end
            end
%-------------遗忘机制：感染节点以delta变为H节点---------------------
            if  s(i,t)==2
                p1=rand(1,1);
                p2=rand(1,1);
%                 if p1<=delta;
%                    s(i,t+1)=3;
%                 elseif p2<=epsilon;
%                     s(i,t+1)=4;
                        if p1<=epsilon;
                   s(i,t+1)=4;
                elseif p2<=delta;
                    s(i,t+1)=3;
                else s(i,t+1)=2;
                end
            end
%-------------被动记忆机制：H节点以ksai变为I节点---------------------
            if  s(i,t)==3
                k=0;                       %统计遗忘节点的邻居节点中是感染节点的数量
                Adjac1=find( A0(i,:)==1 ); %找出节点i的邻居节点位置，就在这里邻接矩阵换成免疫后的新矩阵A0
                num_Adjac1=length(Adjac1); %计算节点i的度
                for g=1:num_Adjac1
                   if s(Adjac1(g),t)==2  %%%%
                      k=k+1;             %统计健康节点的邻居节点中是感染节点的数量k
                   end
                end
                p1=rand(1,1);
                 p2=rand(1,1);
                  e_p=1-(1-xi)^k;   
%                 if p1<=e_p||p2<=eta;
%                    s(i,t+1)=2;
if p1<=e_p
                    s(i,t+1)=2;
                elseif p2<=eta;
                s(i,t+1)=2;
                else s(i,t+1)=3;
                end
            end
%--------------感染节点I变为免疫节点R------------------------
            if s(i,t)==4
                s(i,t+1)=4;
            end
        end
 end
%6.统计谣言每次扩散后，网络中的感染人数与健康人数
 for k2=1:T     
        number_S(z,k2)=length(find(s(:,k2)==0));  %健康节点的密度
        number_E(z,k2)=length(find(s(:,k2)==1));  %感染节点的密度
        number_I(z,k2)=length(find(s(:,k2)==2));  %感染节点的密度
         number_H(z,k2)=length(find(s(:,k2)==3));  %感染节点的密度
        number_R(z,k2)=length(find(s(:,k2)==4));  %免疫节点的密度
    end
    end       
        number_S(z,t)=length(find(s(:,t)==0));              
        number_E(z,t)=length(find(s(:,t)==1));                  
        number_I(z,t)=length(find(s(:,t)==2));                   
        number_H(z,t)=length(find(s(:,t)==3));                   
        number_R(z,t)=length(find(s(:,t)==4));                 
%统计各个节点的平均密度
 number_kS=sum(number_S,1)./s_number./N;
 number_kE=sum(number_E,1)./s_number./N;
 number_kI=sum(number_I,1)./s_number./N;
 number_kH=sum(number_H,1)./s_number./N;
 number_kR=sum(number_R,1)./s_number./N; 
number_basic_I(f,:)= number_kI;
end     
 %画出各节点的变化趋势
 t=1:T;  
 figure;
% plot(t,(number_basic_I(1,:))','-*',t,(number_basic_I(2,:))','-o',t,(number_basic_I(3,:))','-s');
 hold on
% axis([1,40,0,0.5]);

plot(t,(number_basic_I(1,:))','-*','color',[77/256 133/256 189/256]);
plot(t,(number_basic_I(2,:))','-o','color',[247/256 144/256 61/256]);
plot(t,(number_basic_I(3,:))','.-','color',[89/256 169/256 90/256]);

 box on;
 grid off;

xlabel('\fontname{Times new roman}t');ylabel('\fontname{Times new roman}I(t)');  
%%%%% 坐标轴解释                    %%%%% 图例
%  legend('\fontname{Times new roman}\lambda1=0','\fontname{Times new roman}\lambda1=0.1','\fontname{Times new roman}\lambda1=0.3','\fontname{Times new roman}\lambda1=0.5','\fontname{Times new roman}\lambda1=0.5')
 legend('\fontname{Times new roman}\lambda_1=0','\fontname{Times new roman}\lambda_1=0.2','\fontname{Times new roman}\lambda_1=0.4')
         
