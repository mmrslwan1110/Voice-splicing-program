% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   功能：音频格式转换（任意格式互转）           % 
% %   作者：Mr-Ma Technology(马健维)             %
% %   时间：2020.04.09                           %
% %   转载请注明出处                             %
% %   https://blog.csdn.net/qq_29225913/article/details/105445028
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% 支持的转换格式  %%%%%%%%%%%
% % |-----------------------------------------------------|
% % |            格式                     |  输入 |  输出  | 
% % |-----------------------------------------------------|
% % |WAVE (.wav)                          |  支持 |  支持  |
% % |OGG (.ogg)                           |  支持 |  支持  |
% % |FLAC (.flac)                         |  支持 |  支持  |
% % |AU (.au)                             |  支持 | 不支持 |
% % |AIFF(.aiff/.aif)                     |  支持 | 不支持 |
% % |AIFC (.aifc)                         |  支持 | 不支持 |
% % |MP3 (.mp3) (注意电脑版本)             |  支持 | 不支持 |
% % |MPEG-4 AAC(.m4a/.mp4)(注意电脑版本)   |  支持 |  支持  |
% % |-----------------------------------------------------|
% 
% % 更多详细信息，请参考 Matlab 相关描述（通道数、比特率、位深）
% % audioread ：https://ww2.mathworks.cn/help/matlab/ref/audioread.html
% % audiowrite ：https://ww2.mathworks.cn/help/matlab/ref/audiowrite.html?ue&s_tid=gn_loc_dropp
% 
% %% 清空内存里的变量  %%%%%%%%%%%
% clear; 
% close all; 
% clc;
% %%  参数定义  %%%%%%%%%%%
% % 音频文件参数
% % source_FS = 16000; %原采样率（该值弃用，从文件读取）
% dest_FS = 48000; % 目标采样率
% bps = 16;        % 位深 -audiowrite 参数（适用于输出 .wav/.flac 文件）
% qlt = 80;       % Quality -audiowrite 参数（适用于输出 .ogg 文件）
% br  = 0;         % BitRate -audiowrite 参数（适用于输出 .m4a/.mp4 文件） 
% % Tfile = 0;       % 声音片段的总时长s (0为不修改) 暂时没用这个参数
% 
% % 音频通道
% stereo = 'stereo';            % 输出音频通道 'mono'：单声道，'stereo'：立体声
% audio_ch = 'left';          % 立体声输入，单声道输出时，选用哪一通道？ 'left'，'right'，'both'
% Left_channel_ratio = 0.5;   % 左声道占比
% 
% % 文件
% inputfile = ('Recording');   % 输入文件名
% outputfile = [inputfile,'.mp3'];    % 输出文件名
% AudioInfo = audioinfo(inputfile)    % 打印输入音频文件的信息（不加分号）
% 
% %播放设置
% isplay = 1; % 0：播放原音频 1：播放生成音频
% 
% % 绘图
% isdraw = 1;                  % 1：打开绘图， 0：关闭绘图， 计算 FFT 需要消耗时间
% t_axis = [-1 1];             % 时域图显示区域(纵坐标)，按实际来改
% f_axis = [0 20000 0 0.002];  % 频域图显示区域，按实际来改
% 
% %%  原始音频文件读取  %%%%%%%%%%%
% [y_input,source_FS]=audioread(inputfile);         % 将 WAV 文件转换成变量
% fs_str = ['(Fs:',num2str(source_FS/1000),'KHz)']; % fs 的字符串，在图上显示文件 fs
% whos y_input                                      % 打印y的信息
% if isplay == 0                                    % 如果选择播放原音频
%     if source_FS < 192000                         %     电脑小于 192K 采样率才能播放
%         sound(y_input,source_FS)                  %     调用可以播放声音的函数 sound()
%     end
% end
% 
%     
% %%  打印输入波形图像  %%%%%%%%%%%
% if isdraw == 1
%     figure('Name','输入文件波形','NumberTitle','off');
%     [N,n] = size(y_input);                  % 计算文件的长度与声道数量
%     if n==1                                 % 单声道
%         y1 = y_input(:,1);                  % 提取声道数据
%         subplot(2,1,1);             
%         plot(y1);                           % 描绘左声道图像
%         axis([0 N t_axis]);                 % 根据实际情况，限制可显示区间
%         xlabel('采样点');ylabel('幅度');
%         title(['单声道原始波形 ',fs_str]);
% 
%         subplot(2,1,2);
%         NFFT = source_FS * floor(N / source_FS);   % FFT 点数取 fs 的倍数
%         Y_fft = abs(fft(y1,NFFT));                 % 左声道 FFT 计算
%         A_FFT1 = Y_fft * 2 / NFFT;                 % 转换为幅值
%         A_FFT1(1) = 0;                             % 忽略直流分量的影响
%         f = (0:(NFFT/2)) * source_FS / NFFT;       % 频率范围,横坐标（半边谱）
%         plot(f ,A_FFT1(1:NFFT/2 +1));              % 描绘 FFT 图像
%         axis(f_axis);                              % 根据实际情况，限制可显示区间
%         xlabel('频率/Hz');ylabel('振幅');
%         title('单声道原始波形频谱')
% 
%     else                                    % 立体声                          
%         y1 = y_input(:,1);                  % 提取左声道数据
%         y2 = y_input(:,2);                  % 提取右声道数据
% 
%         subplot(2,2,1);                     
%         plot(y1);                            % 描绘左声道图像
%         axis([0 N t_axis]);                  % 根据实际情况，限制可显示区间
%         xlabel('采样点');ylabel('幅度');
%         title(['左声道原始波形 ',fs_str]);
% 
%         subplot(2,2,2);
%         plot(y2);                            % 描绘右声道图像
%         axis([0 N t_axis]);                  % 根据实际情况，限制可显示区间
%         xlabel('采样点');ylabel('幅度');
%         title(['右声道原始波形 ',fs_str]);
% 
%         subplot(2,2,3); 
%         NFFT = source_FS * floor(N / source_FS);    % FFT 点数取 fs 的倍数
%         Y_fft1 = abs(fft(y1,NFFT));                 % 左声道 FFT 计算
%         A_FFT1 = Y_fft1 * 2 / NFFT;                 % 转换为幅值
%         A_FFT1(1) = 0;                              % 忽略直流分量的影响
%         f = (0:(NFFT/2)) * source_FS / NFFT;        % 频率范围,横坐标（半边谱）
% 
%         plot(f ,A_FFT1(1:NFFT/2 +1));               % 描绘 FFT 图像
%         axis(f_axis);                               % 根据实际情况，限制可显示区间
%         xlabel('频率/Hz');ylabel('振幅');
%         title('左声道原始波形频谱')
% 
%         subplot(2,2,4);   
%         NFFT = source_FS * floor(N / source_FS);    % FFT 点数取 fs 的倍数
%         Y_fft2 = abs(fft(y2,NFFT));                 % 左声道 FFT 计算
%         A_FFT2 = Y_fft2 * 2 / NFFT;                 % 转换为幅值
%         A_FFT2(1) = 0;                              % 忽略直流分量的影响
%         f = (0:(NFFT/2)) * source_FS / NFFT;        % 频率范围,横坐标（半边谱）
% 
%         plot(f ,A_FFT2(1:NFFT/2 +1));               % 描绘 FFT 图像
%         axis(f_axis);                               % 根据实际情况，限制可显示区间
%         xlabel('频率/Hz');ylabel('振幅');
%         title('右声道原始波形频谱')
%     end    
% end
% %% 输出参数设置  %%%%%%%%%%%
% if strcmp(stereo,'mono')                % 输出选项为单声道
%     if n == 1                           % 如果输入文件是单声道文件
%         y_ratio = y1;                   % 那么输出文件的数据与输入文件的数据一致
%     else                                % 如果输入是立体声文件
%         if strcmp(audio_ch,'left')      %     输出通道选择 左声道
%             y_ratio = y1;               %     将左声道的数据作为输出数据
%         elseif strcmp(audio_ch,'right') %     输出通道选择 右声道
%             y_ratio = y2;               %     将右声道的数据作为输出数据
%         else                            %     输出通道选择其他选项（both）
%             y_ratio = y1*Left_channel_ratio + y2*(1-Left_channel_ratio); % 就将左声道数据与右声道数据，按照设定的比例来混合。
%         end
%     end
% else                                    % 如果输出选项为立体声
%     if n == 1                           % 输入文件是单声道
%         y_ratio(:,1) = y1;              % 那么左右声道都用输入的数据
%         y_ratio(:,2) = y1;
%     else                                % 输入文件是立体声
%         y_ratio = y_input;              % 那么输出文件的数据与输入文件的数据一致
%     end
% end
% 
% %% 生成输出音频并转码  %%%%%%%%%%%
% y_out = resample(y_ratio,dest_FS,source_FS);              % 对即将输出的数据，进行重采样到设定的采样率
% audiowrite(outputfile,y_out,dest_FS,'BitsPerSample',bps); % 将输出数据存储为音频文件。
% %audiowrite(outputfile,y_out,dest_FS,'Quality',qlt); % 将输出数据存储为音频文件。
% 
% %% 输出音频文件读取  %%%%%%%%%%%
% [y_output,dest_FS]=audioread(outputfile);            % 将上一步生成的文件读取出来保存为变量
% fs_str = ['(Fs:',num2str(dest_FS/1000),'KHz)'];      % 采样率字符串
% 
% if isplay == 1                                       % 如果选择播放生成音频
%     if dest_FS < 192000                              %     电脑小于 192K 采样率才能播放
%         sound(y_output,dest_FS)                      %     调用可以播放声音的函数 sound()
%     end
% end
% 
% %%  打印输出波形图像 （跟打印输入波形图像类似，由于后续复用性较低，就不单独编写函数了）  %%%%%%%%%%%
% if isdraw == 1
%     figure('Name','输出文件波形','NumberTitle','off');    % 打开第2个figure
%     [N,n] = size(y_output);                    % 计算文件的长度与声道数量 
%     if n==1                                    % 单声道
%         y1 = y_output(:,1);                    % 提取声道数据
% 
%         subplot(2,1,1);
%         plot(y1);                              % 描绘左声道图像
%         axis([0 N t_axis]);                    % 根据实际情况，限制可显示区间
%         xlabel('采样点');ylabel('幅度');
%         title(['单声道生成波形 ',fs_str]);
% 
%         subplot(2,1,2);
%         NFFT = dest_FS * floor(N / dest_FS);    % FFT 点数取 fs 的倍数
%         Y_fft = abs(fft(y1,NFFT));              % 左声道 FFT 计算
%         A_FFT1 = Y_fft * 2 / NFFT;              % 转换为幅值
%         A_FFT1(1) = 0;                          % 忽略直流分量的影响
%         f = (0:(NFFT/2)) * dest_FS / NFFT;      % 频率范围,横坐标（半边谱）
%         plot(f ,A_FFT1(1:NFFT/2 +1));           % 描绘 FFT 图像
%         axis(f_axis);                           % 根据实际情况，限制可显示区间
%         xlabel('频率/Hz');ylabel('振幅');
%         title('单声道生成波形频谱');
% 
%     else                                        % 立体声
%         y1 = y_output(:,1);                     % 提取左声道数据
%         y2 = y_output(:,2);                     % 提取右声道数据
% 
%         subplot(2,2,1);
%         plot(y1);                               % 描绘左声道图像
%         axis([0 N t_axis]);                     % 根据实际情况，限制可显示区间
%         xlabel('采样点');ylabel('幅度');
%         title(['左声道生成波形 ',fs_str]);
% 
%         subplot(2,2,2);
%         plot(y2);                               % 描绘右声道图像
%         axis([0 N t_axis]);                     % 根据实际情况，限制可显示区间
%         xlabel('采样点');ylabel('幅度');
%         title(['右声道生成波形 ',fs_str]);
% 
%         subplot(2,2,3);
%         NFFT = dest_FS * floor(N / dest_FS);    % FFT 点数取 fs 的倍数
%         Y_fft1 = abs(fft(y1,NFFT));             % 左声道 FFT 计算
%         A_FFT1 = Y_fft1 * 2 / NFFT;             % 转换为幅值
%         A_FFT1(1) = 0;                          % 忽略直流分量的影响
%         f = (0:(NFFT/2)) * dest_FS / NFFT;      % 频率范围,横坐标（半边谱）
% 
%         plot(f ,A_FFT1(1:NFFT/2 +1));           % 描绘 FFT 图像
%         axis(f_axis);                           % 根据实际情况，限制可显示区间
%         xlabel('频率/Hz');ylabel('振幅');
%         title('左声道生成波形频谱')
% 
%         subplot(2,2,4);                         
%         NFFT = dest_FS * floor(N / dest_FS);    % FFT 点数取 fs 的倍数
%         Y_fft2 = abs(fft(y2,NFFT));             % 左声道 FFT 计算
%         A_FFT2 = Y_fft2 * 2 / NFFT;             % 转换为幅值
%         A_FFT2(1) = 0;                          % 忽略直流分量的影响
%         f = (0:(NFFT/2)) * dest_FS / NFFT;      % 频率范围,横坐标（半边谱）
% 
%         plot(f ,A_FFT2(1:NFFT/2 +1));           % 描绘 FFT 图像
%         axis(f_axis);                           % 根据实际情况，限制可显示区间
%         xlabel('频率/Hz');ylabel('振幅');
%         title('右声道生成波形频谱')
%     end    
% end


clear all;
cd = 'Sound recordings';
waveFiles = dir(fullfile(cd,'*.m4a'));
len = size(waveFiles,1);

Z = [];
for i = 1:len
    fileName = [cd '/' waveFiles(i).name];
    disp(fileName);
    [X, fs] = audioread(fileName);
%     X=X( ( 1 : int32(fs*0.6) ), : );
    
    Z = [Z; X];
end


audiowrite('final.m4a',Z,fs);
