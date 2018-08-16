function varargout = huffman_GUI(varargin)
% HUFFMAN_GUI M-file for huffman_GUI.fig
%      HUFFMAN_GUI, by itself, creates a new HUFFMAN_GUI or raises the existing
%      singleton*.
%
%      H = HUFFMAN_GUI returns the handle to a new HUFFMAN_GUI or the handle to
%      the existing singleton*.
%
%      HUFFMAN_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HUFFMAN_GUI.M with the given input arguments.
%
%      HUFFMAN_GUI('Property','Value',...) creates a new HUFFMAN_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before huffman_GUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to huffman_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help huffman_GUI

% Last Modified by GUIDE v2.5 10-Jan-2014 18:20:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @huffman_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @huffman_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before huffman_GUI is made visible.
function huffman_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to huffman_GUI (see VARARGIN)

% Choose default command line output for huffman_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes huffman_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = huffman_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function xylj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xylj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function xylj_Callback(hObject, eventdata, handles)
% hObject    handle to xylj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xylj as text
%        str2double(get(hObject,'String')) returns contents of xylj as a double


% --- Executes on button press in dqxy.
function dqxy_Callback(hObject, eventdata, handles)
% hObject    handle to dqxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Filename Pathname]=uigetfile({'*.txt'},'请选择文件');       %选择文件框函数
str=[Pathname Filename];      %得到路径和文件名
xinyuan=textread(str,'%s');   %以字符串的形式，读入文本文件
xinyuan=xinyuan{:};            %将cell型数据转换成char
set(handles.xylj,'string',str)   %显示信源的路径
set(handles.xy,'string',xinyuan)  %显示信源


% --- Executes during object creation, after setting all properties.
function xy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function xy_Callback(hObject, eventdata, handles)
% hObject    handle to xy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xy as text
%        str2double(get(hObject,'String')) returns contents of xy as a double


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in bm.
function bm_Callback(hObject, eventdata, handles)
% hObject    handle to bm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;     %清除command window
data=get(handles.xy,'String');     %读入信源
data=uint8(data);         %将信源转化为uint8
[codeword_OK,simbolsout1,fout1,simbolsout2,fout2,zipped,info] = norm2huff(data);  %编码

%排序前
for i=1:length(fout1)
    fout1str{i}=num2str(fout1(i));          %将排序之前的概率转换成cell型数据，方便listbox使用
    str1{i}=char(simbolsout1(i)-1);         %将排序之前的码元转换成cell型数据，方便listbox使用，char将ASCII码转成字符
end
handles.str1=str1;           %得到排序前码元
handles.fout1str=fout1str;   %得到排序前概率

%排序后
for i=1:length(fout2)
    str2{i}=char(simbolsout2(i)-1);      %将排序之后的码元转换成cell型数据，方便listbox使用，char将ASCII码转成字符  
end
handles.str2=str2;           %得到排序后码元               

%得到编码
for i=1:length(codeword_OK)
    codestr{i}=num2str(double(codeword_OK{i}));
    codestr{i}=codestr{i}(find(codestr{i}~=' '));     %去掉字符串里的空格，得到霍夫曼码
end
handles.codestr=codestr;                 %保存霍夫曼码 

handles.zipped=zipped;                  %保存编码的结果
handles.info=info;

guidata(hObject, handles);             %更新handles结构体



% --- Executes on button press in myzl.
function myzl_Callback(hObject, eventdata, handles)
% hObject    handle to myzl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.xsmyzl,'Value',1)
set(handles.xsmyzl,'String',handles.str1);


% --- Executes during object creation, after setting all properties.
function xsmyzl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xsmyzl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in xsmyzl.
function xsmyzl_Callback(hObject, eventdata, handles)
% hObject    handle to xsmyzl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns xsmyzl contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xsmyzl


% --- Executes during object creation, after setting all properties.
function xsglfb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xsglfb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in xsglfb.
function xsglfb_Callback(hObject, eventdata, handles)
% hObject    handle to xsglfb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns xsglfb contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xsglfb


% --- Executes during object creation, after setting all properties.
function xsmypx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xsmypx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in xsmypx.
function xsmypx_Callback(hObject, eventdata, handles)
% hObject    handle to xsmypx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns xsmypx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xsmypx


% --- Executes during object creation, after setting all properties.
function xshfmm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xshfmm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in xshfmm.
function xshfmm_Callback(hObject, eventdata, handles)
% hObject    handle to xshfmm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns xshfmm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xshfmm


% --- Executes during object creation, after setting all properties.
function xsbmjg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xsbmjg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in xsbmjg.
function xsbmjg_Callback(hObject, eventdata, handles)
% hObject    handle to xsbmjg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns xsbmjg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xsbmjg


% --- Executes during object creation, after setting all properties.
function xsym_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xsym (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function xsym_Callback(hObject, eventdata, handles)
% hObject    handle to xsym (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xsym as text
%        str2double(get(hObject,'String')) returns contents of xsym as a double


% --- Executes on button press in ym.
function ym_Callback(hObject, eventdata, handles)
% hObject    handle to ym (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
unzipped=huff2norm(handles.zipped,handles.info);       %译码
yima=char(unzipped);           %将ASCII码转化成字符
set(handles.xsym,'String',yima)          %显示译码结果


% --- Executes on button press in glfb.
function glfb_Callback(hObject, eventdata, handles)
% hObject    handle to glfb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.xsglfb,'Value',1)           %初始化lisbox
set(handles.xsglfb,'String',handles.fout1str);        %在listbox里显示概率分布


% --- Executes on button press in mypx.
function mypx_Callback(hObject, eventdata, handles)
% hObject    handle to mypx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.xsmypx,'Value',1)         %初始化lisbox
set(handles.xsmypx,'String',handles.str2);        %显示排序后的码字


% --- Executes on button press in hfmm.
function hfmm_Callback(hObject, eventdata, handles)
% hObject    handle to hfmm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.xshfmm,'Value',1)            %初始化lisbox
set(handles.xshfmm,'String',handles.codestr)       %lisbox里显示霍夫曼码


% --- Executes on button press in bmjg.
function bmjg_Callback(hObject, eventdata, handles)
% hObject    handle to bmjg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
codestr=handles.codestr;
bianmajieguo=[];
for i=1:length(codestr)
    bianmajieguo=[bianmajieguo codestr{i}];         %将霍夫曼码连成一串
end
set(handles.xsbmjg,'String',bianmajieguo)           %显示编码
    


% --- Executes on button press in qcbm.
function qcbm_Callback(hObject, eventdata, handles)
% hObject    handle to qcbm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%以下是清除编码的工作
set(handles.xylj,'string',[])
set(handles.xy,'string',[])
set(handles.xylj,'string',[])
set(handles.xsmyzl,'string',[])
set(handles.xsglfb,'string',[])
set(handles.xsmypx,'string',[])
set(handles.xshfmm,'string',[])
set(handles.xsbmjg,'string',[])
set(handles.xsym,'string',[])

function f = frequency(vector)
% FREQUENCY   计算元素出现概率
% FREQUENCY(X) 以双精度数组【1*256】的双精度数组返回0~255个元素的概率值

if ~isa(vector,'uint8'),
	error('input argument must be a uint8 vector')
end
%f=repmat(0,1,256);
% 扫描向量
len=length(vector);
for index=0:255; %注意向量的索引参数是从0开始的
    f(index+1)=sum(vector==uint8(index));
end
% 归一化
f=f./len;
%f = histc(vector(:), 0:255); 
%f = f(:)'/sum(f); % always make a row of it

function vector = huff2norm(zipped,info)
%HUFF2NORM   Huffman 解码器
%   For vectors, HUFF2NORM(X,INFO)根据信息结构体info返回向量zipped的解码结果
%  矩阵参数以 X(:)的形式输入
if ~isa(zipped,'uint8'),
	error('input argument must be a uint8 vector')
end

% 产生01序列
len = length(zipped);
string = repmat(uint8(0),1,len.*8);
bitindex = 1:8;
for index = 1:len,
	string(bitindex+8.*(index-1)) = uint8(bitget(zipped(index),bitindex));
end
	
% 调整字符串
string = logical(string(:)');  % make a row of it
len = length(string);
string((len-info.pad+1):end) = [];  % remove 0 padding
len = length(string);

% 解码
weights = 2.^(0:51);
vector = repmat(uint8(0),1,info.length);
vectorindex = 1;
codeindex = 1;
code = 0;
for index = 1:len,
	code = bitset(code,codeindex,string(index));
	codeindex = codeindex+1;
	byte = decode(bitset(code,codeindex),info);
	if byte>0, % 一个码字
		vector(vectorindex) = byte-1;
		codeindex = 1;
		code = 0;
		vectorindex = vectorindex+1;
    end
% % %     fprintf('解码:\n');
% % %     disp(byte);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function byte = decode(code,info)
byte = info.huffcodes(code);

function [codeword_OK,simbolsout1,fout1,simbolsout2,fout2,zipped,info] = norm2huff(vector)
%simbolsout1    排序之前的ASCII码
%fout1          排序之前的对应概率
%simbolsout2    排序之后的ASCII码
%fout2          排序之后的对应概率
%codeword_OK    霍夫曼码


%NORM2HUFF   Huffman 编码器
%  对于输入向量, NORM2HUFF(X) 返回向量的 Huffman编码后的码串
%  矩阵用 X(:) 的形式输入
%  输入限制为uint8格式，输出uint8的序列.
%
%   [...,INFO] 返回解码器需要的结构信息
%   
%      INFO.pad        =最后添加的比特数
%      INFO.huffcodes  = Huffman 码字;
%      INFO.ratio      =压缩率;
%      INFO.length     = 原始数据长度;
%      INFO.maxcodelen = 最大码长;
%

if ~isa(vector,'uint8'),
	error('input argument must be a uint8 vector')
end

% 保证输入为uint8的数据
vector = vector(:)';     %将输入向量转换为行向量
f = frequency(vector);     %计算各元素出现的概率
%f=char(f)
% % % fprintf('各元素出现的概率:\n');
% % % disp(f);
simbols = find(f~=0);  % first value is 1 not 0!!!
f = f(simbols);

fout1=f;    %输出排序前的概率
simbolsout1=simbols;   %输出排序前的码字

% 寻找出现的所有元素
% fprintf('出现的所有元素:\n');
% disp(simbols);
% fprintf('出现的所有元素的概率:\n');
% disp(f)


% 将元素按照出现的概率排列
[f,sortindex] = sort(f);
simbols = simbols(sortindex);
fout2=f;    %输出排序后的概率
simbolsout2=simbols;   %输出排序后的码字
% fprintf('将元素按照出现的概率排列:\n');
% disp(simbols);
% fprintf('对应的概率排列:\n');
% disp(f);

% 产生码字 generate the codewords as the 52 bits of a double
len = length(simbols);
simbols_index = num2cell(1:len);
codeword_tmp = cell(len,1);
while length(f)>1,
	index1 = simbols_index{1};
	index2 = simbols_index{2};
	codeword_tmp(index1) = addnode(codeword_tmp(index1),uint8(0));
	codeword_tmp(index2) = addnode(codeword_tmp(index2),uint8(1));
	f = [sum(f(1:2)) f(3:end)];
  
	simbols_index = [{[index1 index2]} simbols_index(3:end)];
   

   
	% 将数据重新排列，使两个节点的频率尽量与前一个节点的频率相当
   %resort data in order to have the two nodes with lower frequency as first two
	[f,sortindex] = sort(f);
	simbols_index = simbols_index(sortindex);
     
   
end
% 对应相应的元素与码字
codeword = cell(256,1);
codeword(simbols) = codeword_tmp;


% 计算总的字符串长度
len = 0;
for index=1:length(vector),
	len = len+length(codeword{double(vector(index))+1});
end
	
% 产生01序列
string = repmat(uint8(0),1,len);
pointer = 1;
for index=1:length(vector),
	code = codeword{double(vector(index))+1};
	len = length(code);
	string(pointer+(0:len-1)) = code;

	pointer = pointer+len;
      
end

%如果需要，加零
len = length(string);
pad = 8-mod(len,8);
if pad>0,
	string = [string uint8(zeros(1,pad))];
end


% 保存实际有用的码字
codeword = codeword(simbols); 

codeword_OK=codeword;                   %%%%%得到的霍夫曼码



codelen = zeros(size(codeword));
weights = 2.^(0:23);
maxcodelen = 0;
for index = 1:length(codeword),
	len = length(codeword{index});
	if len>maxcodelen,
		maxcodelen = len;
	end
	if len>0,
		code = sum(weights(codeword{index}==1));
		code = bitset(code,len+1);
		codeword{index} = code;
		codelen(index) = len;
	end
end
codeword = [codeword{:}];

   
% 计算压缩 后的向量
cols = length(string)/8;
string = reshape(string,8,cols);
weights = 2.^(0:7);
zipped = uint8(weights*double(string));




% 存储到一个稀疏矩阵
huffcodes = sparse(1,1); % init sparse matrix
for index = 1:numel(codeword),
	huffcodes(codeword(index),1) = simbols(index);
end

% 产生信息结构体
info.pad = pad; 
info.huffcodes = huffcodes;
info.ratio = cols./length(vector);
info.length = length(vector);
info.maxcodelen = maxcodelen;
% fprintf('信息结构体:\n');
% disp(info);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function codeword_new = addnode(codeword_old,item)
codeword_new = cell(size(codeword_old));
for index = 1:length(codeword_old),
	codeword_new{index} = [item codeword_old{index}];
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



