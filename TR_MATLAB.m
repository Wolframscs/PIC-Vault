%% DSC参数辨识
% da1/dt=A1*exp(-E1/R/T)*a1^n1*(1-a1)^m1，da2/dt=A2*exp(-E2/R/T)*a2^n2*(1-a2)^m2
% Q=mx*Hx1*da1/dt+mx*Hx2*da2/dt  , Unit: J/s=W
% Q/mx=Hx1*da1/dt+Hx2*da2/dt     , Unit: W/g    
 
%% START
% clc;clear;close all                                                        % 工作区变量和命令行窗口
% clear global ts tspan y0  T0 R  Hx  sol  P                                 % 清理先前计算的全局变量
% load DSC1.mat 

%% 定义全局变量global 
global ts tspan y0  T0 R  Hx  sol  P dts                                   % 定义全局变量
mx=6e-3;                                                                   % DSC测试中反应物总质量，单位g
Hx=1465;                                                                   % DSC测试中单位质量反应物总释放能量，单位J/g
dTdt=10/60;                                                                % DSC测试中升温速率，单位K/min换算为K/s
Wx=mx*Hx;                                                                  % DSC测试中反应物总释放能量，单位J
T=Z(:,1)+273.15;                                                           % 温度列赋值，且转换为开尔文温标K表示
P=Z(:,2);                                                                  % 产热功率赋值，单位=mW/mg=W/g
T0=T(1);                                                                   % DSC测试起始温度
ts=(T-T0)/10*60;                                                           % 根据温升速率，计算升温时间，单位s
dts=60/10;                                                                 % 时间间隔
tspan=[ts(1) ts(end)];
R=8.314;                                                                   % 通用气体常数赋值
y0=[0 0];                                                                  % y的初始值，1代表未进行反应，0代表反应结束；

%% 参数辨识，拟合 A1 E1 m1，A2 E2 m3 ,Hx1
nvars=8;                                                                   % 待辨识参数个数，此处为A,E,n共3个；
x0=[44.75082122	12.28004124	2.020248685 48.058178	12.45785225	1.799662786 979.4834719  358]; 
x0=x;
% 待辨识参数初始值，从左至右分别为A1,E1,n1,A2,E2,n2,Hx1
% 参数辨识求解器参数的设定
options=optimset("TolFun",1e-12,"MaxFunEvals",1e+20,"TolX",1e-12,"Display","final","PlotFcns","optimplotfval");
% 设置非默认求解器选项及求解，使用fminsearch函数
x=fminsearch(@fitness,x0,options);  
% x=fminimax(@fitness,x0,[],[],[],[],zeros(size(x0)));
% x= ga(@fitness,nvars,[],[],[],[],zeros(nvars,1));
%% 后处理绘图
A1=exp(x(1));E1=exp(x(2)); m1=x(3)                                         % A E 指数放缩，避免参数数量级影响,exp()为较大数；
A2=exp(x(4));E2=exp(x(5)); m2=x(6);                                        % A E 指数放缩，避免参数数量级影响,exp()为较大数；
Hx1=x(7);Hx2=x(8);


a=deval(sol,ts)';                                                          % 计算a1和a2值
da1dt=diff(a(:,1))/6;
da2dt=diff(a(:,2))/6;
da1dt=[0;da1dt];
da2dt=[0;da2dt];
Q_mx=Hx1.*abs(da1dt)+Hx2.*abs(da2dt);                                      % 单位质量产热功率计算，W/g；

figure(2)
subplot(1,3,1)
plot(ts,a(:,1),ts,a(:,2),'--','LineWidth',2);
legend('反应物浓度1：计算','反应物浓度2：计算')
title('Subplot 1: 反应物浓度')
subplot(1,3,2)
plot(ts,da1dt,ts,da2dt,'--','LineWidth',2);
legend('反应物浓度1：计算','反应物浓度2：计算')
title('Subplot 1: 反应物浓度')
subplot(1,3,3)
plot(ts,Q_mx,ts,P,'--','LineWidth',2);
title('Subplot 2: 产热功率')
legend('产热功率：计算','产热功率：测试')
display(x);
Hxcal=sum(dts*Q_mx);
display(Hxcal);
%% 参数辨识拟合
% 参数辨识目标优化函数定义,采用ode45常微分求解器求解
% Q=mx*Hx1*da1/dt+mx*Hx2*da2/dt  , Unit: J/s=W
% Q/mx=Hx1*da1/dt+Hx2*da2/dt,  Unit: W/g    
function ysqrt = fitness(x)                                                % x为参数辨识的向量，定义为A=exp(x(1));E=exp(x(2));m=x(3); n=x(4);
global ts tspan y0  T0 R  Hx  sol  P dts                                   % 定义全局变量
    opts=odeset('RelTol',1e-8,'AbsTol',1e-8);                              % 求解器参数设置,MaxStep可指定最大时间步长，默认为0.1*abs(t0-tf)
    sol=ode45(@(t,y) odefun(t,y,x,R,T0),tspan,y0,opts);                    % t为求解时间范围，0为alpha初始值，opts为求解器设置选项
    a=deval(sol,ts)';                                                      % 计算a1和a2值
    da1dt=diff(a(:,1))/dts;
    da2dt=diff(a(:,2))/dts;
    da1dt=[0;da1dt];
    da2dt=[0;da2dt];
    Hx1=x(7);   Hx2=x(8);                                                     % Hx1 为反应1释放总能量；Hx2为反应2释放总能量；
    Q_mx=Hx1.*abs(da1dt)+Hx2.*abs(da2dt);                                  % 单位质量产热功率计算，W/g；
    ysqrt=sqrt((sum(abs(Q_mx-P).^2))/length(P))                            % 均方根误差计算，fitness函数返回值
end

%% odefun函数定义
% da1/dt=A1*exp(-E1/R/T)*a1^n1*(1-a1)^m1，ode45求解结果为a1
% da2/dt=A2*exp(-E2/R/T)*a2^n2*(1-a2)^m2，ode45求解结果为a2
function dydt=odefun(t,y,x,R,T0)
    dydt=zeros(2,1);
    Ts=T0+t/60*10;                                                         % 根据时间计算当前温度；
    A1=exp(x(1));E1=exp(x(2)); m1=x(3);                                    % A E 指数放缩，避免参数数量级影响,exp()为较大数；
    A2=exp(x(4));E2=exp(x(5)); m2=x(6);                                    % A E 指数放缩，避免参数数量级影响,exp()为较大数；
    dydt(1)=A1*exp(-E1/R./(Ts)).*max(0,(1-y(1))).^m1;
    dydt(2)=A2*exp(-E2/R./(Ts)).*max(0,(1-y(2))).^m2;
end





添加图片注释，不超过 140 字（可选）
