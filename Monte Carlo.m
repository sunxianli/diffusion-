% function [number_aveS,number_aveE,number_aveI,number_aveH,number_aveR]=I_f_r()  %��������Ϊ������ڽӾ����������Ŷ�ǿ�ȣ������̬ʱ��״̬�ڵ�ı���
clc;clear;
%1.���ɳ�ʼ�ڽӾ���
 load('BAnetwork1.mat');  %����WS���缰�ڽӾ���   %�˴��ı䵼������硢�ڽӾ���
A0=A;
% %2.��������
%S������ת��״̬
lambda1=0.4;        %��������S��ΪE�ĸ��� lambda1
lambda2=0.5;       %֧�֡���S��ΪI̬�ĸ���lambda2
lambda3=0.1;       %���ԡ���S��ΪI̬�ĸ��� lambda3
% %�����������ƣ�H��I֮���ת����
delta=0.2;       %��������I��ΪH�ĸ���  delta
xi=0.5;       %�Ӵ����䡪��H������ΪI�ĸ��� ksai
eta=0.5;       %�Է����䡪��H������ΪI�ĸ��� eta
fa=0.8;  %phi
epsilon=0.4;    %epsilon
T=30;                  %TΪҥ������������ɢ��ʱ���ܳ���ʱ�䲽��Ϊ1
s_number=100;          %sΪ������������������ͨ����η�����ƽ��ֵ
%������η�����
number_S=zeros(s_number,T);    
number_E=zeros(s_number,T);
number_I=zeros(s_number,T);
number_H=zeros(s_number,T);
number_R=zeros(s_number,T);
% number_N=zeros(s_number,T);
N=length(A0);        %��ʼ��N
%3.�������
for z=1:s_number
    s=zeros(N,T);       %��ʼ��s��sΪһ����ά���飬���ڴ��ÿ����ɢ�������и����ڵ��״̬��ά��Ҳ��Ҫ��̬�仯
 %���ɳ�ʼ��Ⱦ��---------------------------
ch_q=ceil(rand(1,1)*N);
s(ch_q,1)=2;
%4��Ϣ��ɢ����
 %ÿ����һ����ɢ�󣬶������Ƿ�ﵽ��̬�����ж�
 for t=1:T-1
%      numI=length(find(s(:,t)==2));       
%         if (numI==0)                         %���tʱ�̸�Ⱦ�ڵ�ĸ���Ϊ0���Ͱ�tʱ���Ժ����нڵ�״̬��Ϊ��tʱ����ͬ���������ﵽ����̬
%             for h=t+1:T
%                 s(:,h)=s(:,t);
%             end
%             break;
%         end
        for i=1:N
 %-------------�����ڵ���alpha����Ⱦ��lambda1-3������ת��״̬---------------------
             if  s(i,t)==0
                 k=0;                       %ͳ�ƽ����ڵ���ھӽڵ����Ǹ�Ⱦ�ڵ������
                 Adjac1=find( A0(i,:)==1 ); %�ҳ��ڵ�i���ھӽڵ�λ�ã����������ڽӾ��󻻳����ߺ���¾���A0
                 num_Adjac1=length(Adjac1); %����ڵ�i�Ķ�
                 for g=1:num_Adjac1
                     if s(Adjac1(g),t)==2    %%%%
                        k=k+1;             %ͳ�ƽ����ڵ���ھӽڵ����Ǹ�Ⱦ�ڵ������k
                     end
                 end
                p1=rand(1,1);
                p2=rand(1,1);
                p3=rand(1,1);
                s1=1-(1-lambda1)^k;
                s2=1-(1-lambda2)^k;
                s3=1-(1-lambda3)^k;
                if p1<=s1
                    s(i,t+1)=1;  %S��ΪE
                elseif p2<=s2
                    s(i,t+1)=2;  %S��ΪI̬
                elseif p3<=s3;
                    s(i,t+1)=4;  %S��ΪR̬
                end
             end
%----------------��ԥ����-E��phiֱ�ӱ�ΪI̬(����״̬������)--------------------
            if  s(i,t)==1
                p=rand(1,1);
                if p<=fa
                    s(i,t+1)=2;
                else s(i,t+1)=1;
                end
            end
%-------------�������ƣ���Ⱦ�ڵ���delta��ΪH�ڵ�---------------------
            if  s(i,t)==2
                p1=rand(1,1);
                p2=rand(1,1);
                if p1<=epsilon;
                   s(i,t+1)=4;
                elseif p2<=delta;;
                    s(i,t+1)=3;
                else s(i,t+1)=2;
                end
            end
%-------------����������ƣ�H�ڵ���ksai��ΪI�ڵ�---------------------
            if  s(i,t)==3
                k=0;                       %ͳ�������ڵ���ھӽڵ����Ǹ�Ⱦ�ڵ������
                Adjac1=find( A0(i,:)==1 ); %�ҳ��ڵ�i���ھӽڵ�λ�ã����������ڽӾ��󻻳����ߺ���¾���A0
                num_Adjac1=length(Adjac1); %����ڵ�i�Ķ�
                for g=1:num_Adjac1
                   if s(Adjac1(g),t)==2  %%%%
                      k=k+1;             %ͳ�ƽ����ڵ���ھӽڵ����Ǹ�Ⱦ�ڵ������k
                   end
                end
                p1=rand(1,1);
                 p2=rand(1,1);
                  e_p=1-(1-xi)^k;   
%                 if p1<=e_p||p2<=eta;
              if p1<=e_p
                    s(i,t+1)=2;
                elseif p2<=eta;
                s(i,t+1)=2;
                else s(i,t+1)=3;
                end
            end
%--------------��Ⱦ�ڵ�I��Ϊ���߽ڵ�R------------------------
            if s(i,t)==4
                s(i,t+1)=4;
            end
        end
 end
%6.ͳ��ҥ��ÿ����ɢ�������еĸ�Ⱦ�����뽡������
 for k2=1:T     
        number_S(z,k2)=length(find(s(:,k2)==0));  %�����ڵ���ܶ�
        number_E(z,k2)=length(find(s(:,k2)==1));  %��Ⱦ�ڵ���ܶ�
        number_I(z,k2)=length(find(s(:,k2)==2));  %��Ⱦ�ڵ���ܶ�
         number_H(z,k2)=length(find(s(:,k2)==3));  %��Ⱦ�ڵ���ܶ�
        number_R(z,k2)=length(find(s(:,k2)==4));  %���߽ڵ���ܶ�
    end
    end       
        number_S(z,t)=length(find(s(:,t)==0));              
        number_E(z,t)=length(find(s(:,t)==1));                  
        number_I(z,t)=length(find(s(:,t)==2));                   
        number_H(z,t)=length(find(s(:,t)==3));                   
        number_R(z,t)=length(find(s(:,t)==4));                 
%ͳ�Ƹ����ڵ��ƽ���ܶ�
 number_kS=sum(number_S,1)./s_number./N;
 number_kE=sum(number_E,1)./s_number./N;
 number_kI=sum(number_I,1)./s_number./N;
 number_kH=sum(number_H,1)./s_number./N;
 number_kR=sum(number_R,1)./s_number./N; 

       
 %�������ڵ�ı仯����
 t=1:T;  
 figure;
 hold on
plot(t,number_kS,'-*');
plot(t,number_kE,'-+');
plot(t,number_kI,'-o');
plot(t,number_kH,'.-');
plot(t,number_kR,'-^');
 box on;
 axis([1,30,0,1]);
 grid off;
xlabel('\fontname{Times new roman}t');ylabel('\fontname{Times new roman}Densities of differenent individuals');  
%%%%% ���������                    %%%%% ͼ��
legend('\fontname{Times new roman}S(t)','\fontname{Times new roman}E(t)','\fontname{Times new roman}I(t)','\fontname{Times new roman}H(t)','\fontname{Times new roman}R(t)')
      
