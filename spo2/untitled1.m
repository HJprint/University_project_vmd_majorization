
function varargout = untitled1(varargin)
% UNTITLED1 MATLAB code for untitled1.fig
%      UNTITLED1, by itself, creates a new UNTITLED1 or raises the existing
%      singleton*.
%
%      H = UNTITLED1 returns the handle to a new UNTITLED1 or the handle to
%      the existing singleton*.
%
%      UNTITLED1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED1.M with the given input arguments.
%
%      UNTITLED1('Property','Value',...) creates a new UNTITLED1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled1

% Last Modified by GUIDE v2.5 28-May-2022 13:02:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled1_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled1_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before untitled1 is made visible.
function untitled1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled1 (see VARARGIN)

% Choose default command line output for untitled1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes untitled1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = untitled1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%全局变量
global T;
global fs;



% --- Executes on button press in wt_button.
function wt_button_Callback(hObject, eventdata, handles)
%% 清空环境
global T;
global fs;

%% 参数设置
[signal1,fs]=audioread(T);%加载数据
signal=signal1(:,1);
%signal=load('moni_noise.dat');%加载数据
w=0.9;%权值 将影响PSO 的全局与局部搜优能力， 值较大，全局搜优能力强，局部搜优能力弱;反之，则局部搜优能力增强，而全局搜优能力减弱。
c1=0.1;%加速度，影响收敛速度
c2=0.1;
dim=2;%6维
swarmsize=10;%粒子群规模，表示有10个解的空间
maxiter=10;%最大循环次数，影响时间
minfit=0.001;%最小适应值
vmax=0.01;
vmin=-0.01;
ub=[10,5000];%解向量的最大限制
lb=[2,100];%解向量的最小限制
tau=0;%两处vm d函数用到此参数
DC=0;%两处vmd函数用到此参数
init=1;%两处vmd函数用到此参数
tol=1e-5;%是1x10的-5次方,两处vmd函数用到此参数

%% 种群初始化
swarm=initialization(swarmsize,dim,ub,lb);
vstep=rand(swarmsize,dim)*(vmax-vmin)+vmin;%粒子群速度矩阵
fswarm=zeros(swarmsize,1);%预设空矩阵，存放适应值


%% 个体极值和群体极值
[bestf,bestindex]=min(fswarm);%求得适应值中的最小适应值，和，其所在的序列
gbest=swarm;%暂时的个体最优解为自己
fgbest=fswarm;%暂时的个体最优适应值
zbest=swarm(bestindex,:);%所在序列的对应的解矩阵序列，全局最佳解
fzbest=bestf;%全局最优适应值


%% 迭代寻优
iter=0;
yfitness=zeros(1,maxiter);%1行10列矩阵，存放10个最优值的空间矩阵

while((iter<maxiter)&&(fzbest>minfit))%迭代10次
    for j=1:swarmsize%循环10个粒子
        % 速度更新
        vstep(j,:)=w*vstep(j,:)+c1*rand*(gbest(j,:)-swarm(j,:))+c2*rand*(zbest-swarm(j,:));
        if vstep(j,:)>vmax  
            vstep(j,:)=vmax;%速度限制
        end
        if vstep(j,:)<vmin
            vstep(j,:)=vmin;
        end
        % 位置更新
        swarm(j,:)=swarm(j,:)+vstep(j,:);
        for k=1:dim
            if swarm(j,k)>ub(k)
                swarm(j,k)=ub(k);%位置限制
            end
            if swarm(j,k)<lb(k)
                swarm(j,k)=lb(k);
            end
        end
       
        
        % 返回超出搜索空间边界的搜索代理
        Flag4ub=swarm(j,:)>ub;%两个矩阵比较，是每一个对应位置的值比较
        Flag4lb=swarm(j,:)<lb;
        swarm(j,:)=(swarm(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
         % 适应值
         for ii=1:swarm(j,1)%size(Positions,1)=10，每一个i时，循环Positions第一个参数的x次
            bao=hilbert(u(ii,:));
            bao=abs(bao);
            p=bao./sum(bao);
            e110(ii)=-sum(p.*log10(p))
         end
       fswarm(j,:)=min(e110);%计算每一个position的得分
       %x=swarm(j,:); 
       %fswarm(j,:)=objfun(x,signal);%计算每一个position的得分
        % 可在此处增加约束条件，若满足约束条件，则进行适应值计算
        
        %
        % 个体最优更新
        if fswarm(j)<fgbest(j) %如果当前的函数值比个体最优值小
            gbest(j,:)=swarm(j,:);%个体最优解更新
            fgbest(j)=fswarm(j);%个体最优值更新
        end
        % 群体最优更新
        if fswarm(j)<fzbest%如果当前的函数值比群体最优值大
            zbest=swarm(j,:);%群体最优解更新
            fzbest=fswarm(j);%群体最优值更新
        end
    end
    iter=iter+1;
    yfitness(iter)=fzbest;

end


%% 画图
%% 最优参数
xxx=zbest;
alpha=round(xxx(2));%需要优化的参数
K=round(xxx(1));%需要优化的参数
%% VMD分解
[u, u_hat, omega] = VMD(signal, alpha, tau, K, DC, init, tol);
omega;

%figure(1);%建立幕布
%解决中文字体显示问题
set(0,'defaultAxesFontName', 'Monospaced');
set(0,'defaultAxesFontSize', 10);

plot(handles.wt_tu4,u(1,:),'k'); 
plot(handles.wt_tu5,u(2,:),'k'); 
plot(handles.wt_tu6,u(3,:),'k'); 
plot(handles.wt_tu7,u(4,:),'k'); 
%plot(handles.wt_tu8,u(5,:),'k'); 

signal2=zeros(1,size(signal,1) );
for k=1:K
    if  max(u(k,:))>0.2
        signal2=signal2+u(k,:);
    end
end

fname='ok.wav';
audiowrite(fname,signal2,fs);
[signal2,fs]=audioread(fname);
sound(signal2,fs);

%原始信号 
plot(handles.wt_tu,signal,'k');
set(handles.wt_tu,'XGrid','on','YGrid','on');
%去噪信号 
plot(handles.wt_tu2,signal2,'k');
set(handles.wt_tu2,'XGrid','on','YGrid','on');
%迭代次数
plot(handles.wt_tu3,yfitness,'k');
set(handles.wt_tu3,'XGrid','on','YGrid','on');


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
global T;
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile('.*','选择信号文件');
fpath=[pathname filename];
T=fpath;


% --- Executes on button press in button1.
function button1_Callback(hObject, eventdata, handles)
[signal,fs]=audioread('voice.m4a');%加载数据
sound(signal,fs);
% hObject    handle to button1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button2.
function button2_Callback(hObject, eventdata, handles)
[signal,fs]=audioread('zao.wav');%加载数据
sound(signal,fs);
% hObject    handle to button2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button3.
function button3_Callback(hObject, eventdata, handles)
[signal,fs]=audioread('he.mp3');%加载数据
sound(signal,fs);
% hObject    handle to button3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button4.
function button4_Callback(hObject, eventdata, handles)
[signal,fs]=audioread('ok.wav');%加载数据
sound(signal,fs);
% hObject    handle to button4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
