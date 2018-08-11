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
[Filename Pathname]=uigetfile({'*.txt'},'��ѡ���ļ�');       %ѡ���ļ�����
str=[Pathname Filename];      %�õ�·�����ļ���
xinyuan=textread(str,'%s');   %���ַ�������ʽ�������ı��ļ�
xinyuan=xinyuan{:};            %��cell������ת����char
set(handles.xylj,'string',str)   %��ʾ��Դ��·��
set(handles.xy,'string',xinyuan)  %��ʾ��Դ


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
clc;     %���command window
data=get(handles.xy,'String');     %������Դ
data=uint8(data);         %����Դת��Ϊuint8
[codeword_OK,simbolsout1,fout1,simbolsout2,fout2,zipped,info] = norm2huff(data);  %����

%����ǰ
for i=1:length(fout1)
    fout1str{i}=num2str(fout1(i));          %������֮ǰ�ĸ���ת����cell�����ݣ�����listboxʹ��
    str1{i}=char(simbolsout1(i)-1);         %������֮ǰ����Ԫת����cell�����ݣ�����listboxʹ�ã�char��ASCII��ת���ַ�
end
handles.str1=str1;           %�õ�����ǰ��Ԫ
handles.fout1str=fout1str;   %�õ�����ǰ����

%�����
for i=1:length(fout2)
    str2{i}=char(simbolsout2(i)-1);      %������֮�����Ԫת����cell�����ݣ�����listboxʹ�ã�char��ASCII��ת���ַ�  
end
handles.str2=str2;           %�õ��������Ԫ               

%�õ�����
for i=1:length(codeword_OK)
    codestr{i}=num2str(double(codeword_OK{i}));
    codestr{i}=codestr{i}(find(codestr{i}~=' '));     %ȥ���ַ�����Ŀո񣬵õ���������
end
handles.codestr=codestr;                 %����������� 

handles.zipped=zipped;                  %�������Ľ��
handles.info=info;

guidata(hObject, handles);             %����handles�ṹ��



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
unzipped=huff2norm(handles.zipped,handles.info);       %����
yima=char(unzipped);           %��ASCII��ת�����ַ�
set(handles.xsym,'String',yima)          %��ʾ������


% --- Executes on button press in glfb.
function glfb_Callback(hObject, eventdata, handles)
% hObject    handle to glfb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.xsglfb,'Value',1)           %��ʼ��lisbox
set(handles.xsglfb,'String',handles.fout1str);        %��listbox����ʾ���ʷֲ�


% --- Executes on button press in mypx.
function mypx_Callback(hObject, eventdata, handles)
% hObject    handle to mypx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.xsmypx,'Value',1)         %��ʼ��lisbox
set(handles.xsmypx,'String',handles.str2);        %��ʾ����������


% --- Executes on button press in hfmm.
function hfmm_Callback(hObject, eventdata, handles)
% hObject    handle to hfmm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.xshfmm,'Value',1)            %��ʼ��lisbox
set(handles.xshfmm,'String',handles.codestr)       %lisbox����ʾ��������


% --- Executes on button press in bmjg.
function bmjg_Callback(hObject, eventdata, handles)
% hObject    handle to bmjg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
codestr=handles.codestr;
bianmajieguo=[];
for i=1:length(codestr)
    bianmajieguo=[bianmajieguo codestr{i}];         %��������������һ��
end
set(handles.xsbmjg,'String',bianmajieguo)           %��ʾ����
    


% --- Executes on button press in qcbm.
function qcbm_Callback(hObject, eventdata, handles)
% hObject    handle to qcbm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%�������������Ĺ���
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
% FREQUENCY   ����Ԫ�س��ָ���
% FREQUENCY(X) ��˫�������顾1*256����˫�������鷵��0~255��Ԫ�صĸ���ֵ

if ~isa(vector,'uint8'),
	error('input argument must be a uint8 vector')
end
%f=repmat(0,1,256);
% ɨ������
len=length(vector);
for index=0:255; %ע�����������������Ǵ�0��ʼ��
    f(index+1)=sum(vector==uint8(index));
end
% ��һ��
f=f./len;
%f = histc(vector(:), 0:255); 
%f = f(:)'/sum(f); % always make a row of it

function vector = huff2norm(zipped,info)
%HUFF2NORM   Huffman ������
%   For vectors, HUFF2NORM(X,INFO)������Ϣ�ṹ��info��������zipped�Ľ�����
%  ��������� X(:)����ʽ����
if ~isa(zipped,'uint8'),
	error('input argument must be a uint8 vector')
end

% ����01����
len = length(zipped);
string = repmat(uint8(0),1,len.*8);
bitindex = 1:8;
for index = 1:len,
	string(bitindex+8.*(index-1)) = uint8(bitget(zipped(index),bitindex));
end
	
% �����ַ���
string = logical(string(:)');  % make a row of it
len = length(string);
string((len-info.pad+1):end) = [];  % remove 0 padding
len = length(string);

% ����
weights = 2.^(0:51);
vector = repmat(uint8(0),1,info.length);
vectorindex = 1;
codeindex = 1;
code = 0;
for index = 1:len,
	code = bitset(code,codeindex,string(index));
	codeindex = codeindex+1;
	byte = decode(bitset(code,codeindex),info);
	if byte>0, % һ������
		vector(vectorindex) = byte-1;
		codeindex = 1;
		code = 0;
		vectorindex = vectorindex+1;
    end
% % %     fprintf('����:\n');
% % %     disp(byte);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function byte = decode(code,info)
byte = info.huffcodes(code);

function [codeword_OK,simbolsout1,fout1,simbolsout2,fout2,zipped,info] = norm2huff(vector)
%simbolsout1    ����֮ǰ��ASCII��
%fout1          ����֮ǰ�Ķ�Ӧ����
%simbolsout2    ����֮���ASCII��
%fout2          ����֮��Ķ�Ӧ����
%codeword_OK    ��������


%NORM2HUFF   Huffman ������
%  ������������, NORM2HUFF(X) ���������� Huffman�������봮
%  ������ X(:) ����ʽ����
%  ��������Ϊuint8��ʽ�����uint8������.
%
%   [...,INFO] ���ؽ�������Ҫ�Ľṹ��Ϣ
%   
%      INFO.pad        =�����ӵı�����
%      INFO.huffcodes  = Huffman ����;
%      INFO.ratio      =ѹ����;
%      INFO.length     = ԭʼ���ݳ���;
%      INFO.maxcodelen = ����볤;
%

if ~isa(vector,'uint8'),
	error('input argument must be a uint8 vector')
end

% ��֤����Ϊuint8������
vector = vector(:)';     %����������ת��Ϊ������
f = frequency(vector);     %�����Ԫ�س��ֵĸ���
%f=char(f)
% % % fprintf('��Ԫ�س��ֵĸ���:\n');
% % % disp(f);
simbols = find(f~=0);  % first value is 1 not 0!!!
f = f(simbols);

fout1=f;    %�������ǰ�ĸ���
simbolsout1=simbols;   %�������ǰ������

% Ѱ�ҳ��ֵ�����Ԫ��
% fprintf('���ֵ�����Ԫ��:\n');
% disp(simbols);
% fprintf('���ֵ�����Ԫ�صĸ���:\n');
% disp(f)


% ��Ԫ�ذ��ճ��ֵĸ�������
[f,sortindex] = sort(f);
simbols = simbols(sortindex);
fout2=f;    %��������ĸ���
simbolsout2=simbols;   %�������������
% fprintf('��Ԫ�ذ��ճ��ֵĸ�������:\n');
% disp(simbols);
% fprintf('��Ӧ�ĸ�������:\n');
% disp(f);

% �������� generate the codewords as the 52 bits of a double
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
   

   
	% �������������У�ʹ�����ڵ��Ƶ�ʾ�����ǰһ���ڵ��Ƶ���൱
   %resort data in order to have the two nodes with lower frequency as first two
	[f,sortindex] = sort(f);
	simbols_index = simbols_index(sortindex);
     
   
end
% ��Ӧ��Ӧ��Ԫ��������
codeword = cell(256,1);
codeword(simbols) = codeword_tmp;


% �����ܵ��ַ�������
len = 0;
for index=1:length(vector),
	len = len+length(codeword{double(vector(index))+1});
end
	
% ����01����
string = repmat(uint8(0),1,len);
pointer = 1;
for index=1:length(vector),
	code = codeword{double(vector(index))+1};
	len = length(code);
	string(pointer+(0:len-1)) = code;

	pointer = pointer+len;
      
end

%�����Ҫ������
len = length(string);
pad = 8-mod(len,8);
if pad>0,
	string = [string uint8(zeros(1,pad))];
end


% ����ʵ�����õ�����
codeword = codeword(simbols); 

codeword_OK=codeword;                   %%%%%�õ��Ļ�������



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

   
% ����ѹ�� �������
cols = length(string)/8;
string = reshape(string,8,cols);
weights = 2.^(0:7);
zipped = uint8(weights*double(string));




% �洢��һ��ϡ�����
huffcodes = sparse(1,1); % init sparse matrix
for index = 1:numel(codeword),
	huffcodes(codeword(index),1) = simbols(index);
end

% ������Ϣ�ṹ��
info.pad = pad; 
info.huffcodes = huffcodes;
info.ratio = cols./length(vector);
info.length = length(vector);
info.maxcodelen = maxcodelen;
% fprintf('��Ϣ�ṹ��:\n');
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



