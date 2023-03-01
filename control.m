clc,clear;
% 参数设置
test = -1;%测试变量，循环条件
delta1 = 0.001;%用于比较最后的精度
N = 2000;%节点的个数
T=15;%时间长度
t=linspace(0,T,N+1);
h=T/N;
h2 = h/2;
%定义传播参数
lambda1=0.4;lambda2=0.5;lambda3=0.1;
epsilon=0.4;delta=0.2;
eta=0.5;xi=0.5;fa=0.8;
 u_max=1;v_max=1;w_max=1;
Du = [6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	35	36	37	38	39	40	41	42	44	45	46	48	50	51	52	53	54	57	58	59	61	62	63	67	70	71	77	78	84	86	92	94	96	99	120	126	137	147	152	186
]; % 度
Pk = [0.243	0.16	0.104	0.0905	0.0695	0.0555	0.04	0.031	0.0235	0.018	0.019	0.0165	0.014	0.0085	0.0125	0.0095	0.0035	0.004	0.006	0.0045	0.0055	0.0055	0.005	0.0035	0.0035	0.002	0.001	0.0025	0.002	0.002	0.0035	0.0015	0.001	0.0015	0.001	0.002	0.001	0.002	0.001	0.0015	0.001	0.0015	0.001	0.001	0.0005	0.0005	0.001	0.0005	0.001	0.0005	0.0005	0.0005	0.0005	0.0005	0.0005	0.0005	0.0005	0.001	0.001	0.0005	0.0005	0.0005	0.0005	0.0005	0.0005	0.0005	0.0005	0.0005
];
 A=0.05*ones(length(Du),1);
 B=0.8*ones(length(Du),1);
 C=0.5*ones(length(Du),1);
%平均度计算
k2 = 0.0;
for i=1:length(Du)
   k2 = k2+Du(i)*Pk(i); 
end
% 控制变量的初始猜测值
u = zeros(N+1,length(Du));        %创建N+1行Du_Num列的0矩阵
v = zeros(N+1,length(Du));
w = zeros(N+1,length(Du));
% 为状态变量分配内存
S = zeros(N+1,length(Du));         %创建M+1行Du_Num列的0矩阵
E = zeros(N+1,length(Du));
I = zeros(N+1,length(Du));
H = zeros(N+1,length(Du));
R = zeros(N+1,length(Du));
lambdaS = zeros(N+1,length(Du));
lambdaE = zeros(N+1,length(Du));
lambdaI = zeros(N+1,length(Du));
lambdaH = zeros(N+1,length(Du));
lambdaR = zeros(N+1,length(Du));
% 为状态变量赋初值
S(1,:) = 0.8*ones(length(Du),1);     %S矩阵的第一行赋值
E(1,:) = 0.05*ones(length(Du),1);
I(1,:) = 0.1*ones(length(Du),1);
H(1,:) = 0.05*ones(length(Du),1);
R(1,:) = 0*ones(length(Du),1);
theta = zeros(N+1,1);%定义与i的接触矩阵
while (test < 0 )
    oldu = u;
    oldv = v;
    oldw = w;
    oldS = S;
    oldE =E;
    oldI = I;
    oldH = H;
    oldR = R;
    oldlambdaS = lambdaS;
    oldlambdaE = lambdaE;
    oldlambdaI = lambdaI;
    oldlambdaH = lambdaH;
    oldlambdaR = lambdaR;
    %前向求解约束方程
    for t=1:N        %每一步
      %计算与i接触的概率和  
       theta(t) = 0.0;
       for i=1:length(Du)
           theta(t) = theta(t)+Du(i)*Pk(i)*I(t,i);
       end
       theta(t) = theta(t)/k2;
       
       for k=1:length(Du)     %求解每个度的状态
            %四阶龙格库塔法
        m11 =-Du(k)*S(t,k)*theta(t);
        m12 = lambda1*Du(k)*S(t,k)*theta(t)-fa*E(t,k)-u(t,k)*E(t,k);
        m13 = lambda2*Du(k)*S(t,k)*theta(t)+fa*E(t,k)+eta*H(t,k)+xi*Du(k)*H(t,k)*theta(t)-delta*I(t,k)-epsilon*I(t,k)-v(t,k)*I(t,k);
        m14 = delta*I(t,k)-eta*H(t,k)-xi*H(t,k)*Du(k)*theta(t)-w(t,k)*H(t,k);
        m15=lambda3*Du(k)*S(t,k)*theta(t)+epsilon*I(t,k)+u(t,k)*E(t,k)+v(t,k)*I(t,k)+w(t,k)*H(t,k);

        m21 =-Du(k)*(S(t,k)+h2*m11)*theta(t);
        m22 = lambda1*Du(k)*(S(t,k)+h2*m11)*theta(t)-fa*(E(t,k)+h2*m12)-0.5*(u(t,k)+u(t+1,k))*(E(t,k)+h2*m12);
        m23 = lambda2*Du(k)*(S(t,k)+h2*m11)*theta(t)+fa*(E(t,k)+h2*m12)+eta*(H(t,k)+h2*m14)+xi*Du(k)*(H(t,k)+h2*m14)*theta(t)-delta*(I(t,k)+h2*m13)-epsilon*(I(t,k)+h2*m13)-0.5*(v(t,k)+v(t+1,k))*(I(t,k)+h2*m13);
        m24 = delta*(I(t,k)+h2*m13)-eta*(H(t,k)+h2*m14)-xi*(H(t,k)+h2*m14)*Du(k)*theta(t)-0.5*(w(t,k)+w(t+1,k))*(H(t,k)+h2*m14);
        m25=lambda3*Du(k)*(S(t,k)+h2*m11)*theta(t)+epsilon*(I(t,k)+h2*m13)+0.5*(u(t,k)+u(t+1,k))*(E(t,k)+h2*m12)+0.5*(v(t,k)+v(t+1,k))*(I(t,k)+h2*m13)+0.5*(w(t,k)+w(t+1,k))*(H(t,k)+h2*m14);   
           
        m31 =-Du(k)*(S(t,k)+h2*m21)*theta(t);
        m32 = lambda1*Du(k)*(S(t,k)+h2*m21)*theta(t)-fa*(E(t,k)+h2*m22)-0.5*(u(t,k)+u(t+1,k))*(E(t,k)+h2*m22);
        m33 = lambda2*Du(k)*(S(t,k)+h2*m21)*theta(t)+fa*(E(t,k)+h2*m22)+eta*(H(t,k)+h2*m24)+xi*Du(k)*(H(t,k)+h2*m24)*theta(t)-delta*(I(t,k)+h2*m23)-epsilon*(I(t,k)+h2*m23)-0.5*(v(t,k)+v(t+1,k))*(I(t,k)+h2*m23);
        m34 = delta*(I(t,k)+h2*m23)-eta*(H(t,k)+h2*m24)-xi*(H(t,k)+h2*m24)*Du(k)*theta(t)-0.5*(w(t,k)+w(t+1,k))*(H(t,k)+h2*m24);
        m35=lambda3*Du(k)*(S(t,k)+h2*m21)*theta(t)+epsilon*(I(t,k)+h2*m23)+0.5*(u(t,k)+u(t+1,k))*(E(t,k)+h2*m22)+0.5*(v(t,k)+v(t+1,k))*(I(t,k)+h2*m23)+0.5*(w(t,k)+w(t+1,k))*(H(t,k)+h2*m24);   
           
         
        m41 =-Du(k)*(S(t,k)+h*m31)*theta(t);
        m42 = lambda1*Du(k)*(S(t,k)+h*m31)*theta(t)-fa*(E(t,k)+h*m32)-u(t+1,k)*(E(t,k)+h*m32);
        m43 = lambda2*Du(k)*(S(t,k)+h*m31)*theta(t)+fa*(E(t,k)+h*m32)+eta*(H(t,k)+h*m34)+xi*Du(k)*(H(t,k)+h*m34)*theta(t)-delta*(I(t,k)+h*m33)-epsilon*(I(t,k)+h*m33)-v(t+1,k)*(I(t,k)+h*m33);
        m44 = delta*(I(t,k)+h*m33)-eta*(H(t,k)+h*m34)-xi*(H(t,k)+h*m34)*Du(k)*theta(t)-w(t+1,k)*(H(t,k)+h*m34);
        m45=lambda3*Du(k)*(S(t,k)+h*m31)*theta(t)+epsilon*(I(t,k)+h*m33)+u(t+1,k)*(E(t,k)+h*m32)+v(t+1,k)*(I(t,k)+h*m33)+w(t+1,k)*(H(t,k)+h*m34);   
           
         
            S(t+1,k) = S(t,k)+h/6*(m11 + 2*m21 + 2*m31 + m41);
            E(t+1,k) = E(t,k)+h/6*(m12 + 2*m22 + 2*m32 + m42);
            I(t+1,k) = I(t,k)+h/6*(m13 + 2*m23 + 2*m33 + m43);
            H(t+1,k) =H(t,k)+h/6*(m14 + 2*m24 + 2*m34 + m44);
           R(t+1,k) = R(t,k)+h/6*(m15 + 2*m25 + 2*m35 + m45);
           end
    end
    
    %向后更新
    for t=N+1:-1:2         %每一步
       for k=1:length(Du)     %求解每个度的状态
           theta(t) = 0.0;
           for i=1:length(Du)
               theta(t) = theta(t)+Du(i)*Pk(i)*I(t,i);
           end
           theta(t) = theta(t)/k2;
           %lamda方程
        m11 = Du(k)*lambdaS(t,k)*theta(t)-Du(k)*lambdaE(t,k)*lambda1*theta(t)-Du(k)*lambdaI(t,k)*lambda2*theta(t)-Du(k)*lambdaR(t,k)*lambda3*theta(t);
        m12 = lambdaE(t,k)*(fa+u(t,k))-lambdaI(t,k)*fa-lambdaR(t,k)*u(t,k);
        m13 =-1+ lambdaS(t,k)*S(t,k)*Du(k)^2*Pk(k)/k2-lambdaE(t,k)*lambda1*S(t,k)*Du(k)^2*Pk(k)/k2-lambdaI(t,k)*(lambda2*S(t,k)*Du(k)^2*Pk(k)/k2+xi*H(t,k)*Du(k)^2*Pk(k)/k2-delta-epsilon-v(t,k))-lambdaH(t,k)*(delta-xi*H(t,k)*Du(k)^2*Pk(k)/k2)-lambdaR(t,k)*(lambda3*S(t,k)*Du(k)^2*Pk(k)/k2+epsilon+v(t,k));
        m14 =-lambdaI(t,k)*(eta+xi*Du(k)*theta(t))+lambdaH(t,k)*(w(t,k)+xi*Du(k)*theta(t)+eta)-lambdaR(t,k)*w(t,k);
        m15 =0;
           
        m21 = Du(k)*(lambdaS(t,k)-h2*m11)*theta(t)-Du(k)*(lambdaE(t,k)-h2*m12)*lambda1*theta(t)-Du(k)*(lambdaI(t,k)-h2*m13)*lambda2*theta(t)-Du(k)*(lambdaR(t,k)-h2*m15)*lambda3*theta(t);
        m22 = (lambdaE(t,k)-h2*m12)*(fa+0.5*(u(t,k)+u(t-1,k)))-(lambdaI(t,k)-h2*m13)*fa-(lambdaR(t,k)-h2*m15)*0.5*(u(t,k)+u(t-1,k));
        m23 =-1+ (lambdaS(t,k)-h2*m11)*0.5*(S(t,k)+S(t-1,k))*Du(k)^2*Pk(k)/k2-(lambdaE(t,k)-h2*m12)*lambda1*0.5*(S(t,k)+S(t-1,k))*Du(k)^2*Pk(k)/k2-(lambdaI(t,k)-h2*m13)*(lambda2*0.5*(S(t,k)+S(t-1,k))*Du(k)^2*Pk(k)/k2+xi*0.5*(H(t,k)+H(t-1,k))*Du(k)^2*Pk(k)/k2-delta-epsilon-0.5*(v(t,k)+v(t-1,k)))-(lambdaH(t,k)-h2*m14)*(delta-xi*0.5*(H(t,k)+H(t-1,k))*Du(k)^2*Pk(k)/k2)-(lambdaR(t,k)-h2*m15)*(lambda3*0.5*(S(t,k)+S(t-1,k))*Du(k)^2*Pk(k)/k2+epsilon+0.5*(v(t,k)+v(t-1,k)));
        m24 =-(lambdaI(t,k)-h2*m13)*(eta+xi*Du(k)*theta(t))+(lambdaH(t,k)-h2*m14)*(0.5*(w(t,k)+w(t-1,k))+xi*Du(k)*theta(t)+eta)-(lambdaR(t,k)-h2*m15)*0.5*(w(t,k)+w(t-1,k));
        m25 =0;
           
        m31 = Du(k)*(lambdaS(t,k)-h2*m21)*theta(t)-Du(k)*(lambdaE(t,k)-h2*m22)*lambda1*theta(t)-Du(k)*(lambdaI(t,k)-h2*m23)*lambda2*theta(t)-Du(k)*(lambdaR(t,k)-h2*m25)*lambda3*theta(t);
        m32 = (lambdaE(t,k)-h2*m22)*(fa+0.5*(u(t,k)+u(t-1,k)))-(lambdaI(t,k)-h2*m23)*fa-(lambdaR(t,k)-h2*m25)*0.5*(u(t,k)+u(t-1,k));
        m33 =-1+ (lambdaS(t,k)-h2*m21)*0.5*(S(t,k)+S(t-1,k))*Du(k)^2*Pk(k)/k2-(lambdaE(t,k)-h2*m22)*lambda1*0.5*(S(t,k)+S(t-1,k))*Du(k)^2*Pk(k)/k2-(lambdaI(t,k)-h2*m23)*(lambda2*0.5*(S(t,k)+S(t-1,k))*Du(k)^2*Pk(k)/k2+xi*0.5*(H(t,k)+H(t-1,k))*Du(k)^2*Pk(k)/k2-delta-epsilon-0.5*(v(t,k)+v(t-1,k)))-(lambdaH(t,k)-h2*m24)*(delta-xi*0.5*(H(t,k)+H(t-1,k))*Du(k)^2*Pk(k)/k2)-(lambdaR(t,k)-h2*m25)*(lambda3*0.5*(S(t,k)+S(t-1,k))*Du(k)^2*Pk(k)/k2+epsilon+0.5*(v(t,k)+v(t-1,k)));
        m34 =-(lambdaI(t,k)-h2*m23)*(eta+xi*Du(k)*theta(t))+(lambdaH(t,k)-h2*m24)*(0.5*(w(t,k)+w(t-1,k))+xi*Du(k)*theta(t)+eta)-(lambdaR(t,k)-h2*m25)*0.5*(w(t,k)+w(t-1,k));
        m35 =0;
        
        m41 = Du(k)*(lambdaS(t,k)-h*m31)*theta(t)-Du(k)*(lambdaE(t,k)-h*m32)*lambda1*theta(t)-Du(k)*(lambdaI(t,k)-h*m33)*lambda2*theta(t)-Du(k)*(lambdaR(t,k)-h*m35)*lambda3*theta(t);
        m42 = (lambdaE(t,k)-h*m32)*(fa+u(t-1,k))-(lambdaI(t,k)-h*m33)*fa-(lambdaR(t,k)-h*m35)*u(t-1,k);
        m43 =-1+ (lambdaS(t,k)-h*m31)*S(t-1,k)*Du(k)^2*Pk(k)/k2-(lambdaE(t,k)-h*m32)*lambda1*S(t-1,k)*Du(k)^2*Pk(k)/k2-(lambdaI(t,k)-h*m33)*(lambda2*S(t-1,k)*Du(k)^2*Pk(k)/k2+xi*H(t-1,k)*Du(k)^2*Pk(k)/k2-delta-epsilon-v(t-1,k))-(lambdaH(t,k)-h*m34)*(delta-xi*H(t-1,k)*Du(k)^2*Pk(k)/k2)-(lambdaR(t,k)-h*m35)*(lambda3*S(t-1,k)*Du(k)^2*Pk(k)/k2+epsilon+v(t-1,k));
        m44 =-(lambdaI(t,k)-h*m33)*(eta+xi*Du(k)*theta(t))+(lambdaH(t,k)-h*m34)*(w(t-1,k)+xi*Du(k)*theta(t)+eta)-(lambdaR(t,k)-h*m35)*w(t-1,k);
        m45 =0;
            lambdaS(t-1,k) =lambdaS(t,k) - h/6*(m11 + 2*m21 + 2*m31 + m41);
            lambdaE(t-1,k) = lambdaE(t,k) - h/6*(m12 + 2*m22 + 2*m32 + m42);
            lambdaI(t-1,k) = lambdaI(t,k) - h/6*(m13 + 2*m23 + 2*m33 + m43);
            lambdaH(t-1,k) = lambdaH(t,k) - h/6*(m14 + 2*m24 + 2*m34 + m44);
            lambdaR(t-1,k) = lambdaR(t,k) - h/6*(m15 + 2*m25 + 2*m35 + m45);
            
       end
        
    end 
    
    for t=1:N+1 
        for k=1:length(Du)
            psi1(t,k)=(lambdaE(t,k)-lambdaR(t,k))*E(t,k)/A(k);
            psi2(t,k)=(lambdaI(t,k)-lambdaR(t,k))*I(t,k)/B(k);  
            psi3(t,k)=(lambdaH(t,k)-lambdaR(t,k))*H(t,k)/C(k);
            u1(t,k) = min(max(psi1(t,k),0),u_max);
            v1(t,k) = min(max(psi2(t,k),0),v_max);
            w1(t,k)=min(max(psi3(t,k),0),w_max);
            u(t,k) = 0.5*(u1(t,k) + oldu(t,k)); 
            v(t,k) = 0.5*(v1(t,k) + oldv(t,k)); 
            w(t,k)=0.5*(w1(t,k) + oldw(t,k));  
        end
    end

    % 确定是否循环继续求解最优系统
    temp(1) = delta1*sum(abs(u(:))) - sum(abs(oldu(:) - u(:)));
    temp(2) = delta1*sum(abs(v(:))) - sum(abs(oldv(:) - v(:)));
    temp(3) = delta1*sum(abs(w(:))) - sum(abs(oldw(:) - w(:)));
    temp(4) = delta1*sum(abs(S(:))) - sum(abs(oldS(:) - S(:)));
    temp(5) = delta1*sum(abs(E(:))) - sum(abs(oldE(:) - E(:)));
    temp(6) = delta1*sum(abs(I(:))) - sum(abs(oldI(:) - I(:)));
    temp(7) = delta1*sum(abs(H(:))) - sum(abs(oldH(:) - H(:)));
    temp(8) = delta1*sum(abs(R(:))) - sum(abs(oldR(:) - R(:)));
    temp(9) = delta1*sum(abs(lambdaS(:))) - sum(abs(oldlambdaS(:) - lambdaS(:)));
    temp(10) = delta1*sum(abs(lambdaE(:))) - sum(abs(oldlambdaE(:) - lambdaE(:)));
    temp(11) = delta1*sum(abs(lambdaI(:))) - sum(abs(oldlambdaI(:) - lambdaI(:)));
    temp(12) = delta1*sum(abs(lambdaH(:))) - sum(abs(oldlambdaH(:) - lambdaH(:)));
    temp(13) = delta1*sum(abs(lambdaR(:))) - sum(abs(oldlambdaR(:) - lambdaR(:)));
    test = min(temp);
end
 uaverage=mean(u(:))
 vaverage=mean(v(:))
 waverage=mean(w(:))
figure(1);
J= 0.0;
for k=1:length(Du)
    temp_J = I(:,k)+1/2*A(k)*u(:,k).^2+1/2*B(k)*v(:,k).^2+1/2*C(k)*w(:,k).^2;
    J = J+temp_J;
end
J1=trapz(t,J)
bar(J1)
hold on;
 box on;
 grid off;
%画出五类节点的状态变化图
figure(2);
 for k=1:length(Du);
            S(:,k)= S(:,k)*Pk(k);
            E(:,k)= E(:,k)*Pk(k);
             I(:,k)=I(:,k)*Pk(k);
             H(:,k)= H(:,k)*Pk(k);
             R(:,k)= R(:,k)*Pk(k);
        end
        S=sum(S,2);
        E=sum(E,2);
        I=sum(I,2);
        H=sum(H,2);
        R=sum(R,2);
        
 t=0:0.0075:15;
 hold on
 box on;
 grid off;
% xlim([0,15]);
 plot(t,S','-*',t,E','-+',t,I','-o',t,H','.-',t,R','-^','MarkerIndices',1:80:2001)
xlabel('\fontname{Times new roman}t');ylabel('\fontname{Times new roman}Densities of differenent individuals');  
legend('\fontname{Times new roman}S(t)','\fontname{Times new roman}E(t)','\fontname{Times new roman}I(t)','\fontname{Times new roman}H(t)','\fontname{Times new roman}R(t)')
% figure(3);

% %  plot(t,u(:,1)',t,u(:,20)',t,u(:,40)',t,u(:,68));
%  U=sum(u,2)/length(Du);%行求和   
% W=sum(w,2)/length(Du);%行求和   
%  V=sum(v,2)/length(Du);%行求和 
%  hold on;
%   box on;
%  grid off;
%  t=0:0.0075:15;
% plot(t,U','-*','MarkerIndices',1:80:2001);
% plot(t,V','-O','MarkerIndices',1:80:2001);
% plot(t,W','.-','MarkerIndices',1:80:2001);
% xlabel('\fontname{Times new roman}t');ylabel('\fontname{Times new roman}control strength ');  
% legend('\fontname{Times new roman}<u(t)>','\fontname{Times new roman}<v(t)>','\fontname{Times new roman}<w(t)>')
% 
