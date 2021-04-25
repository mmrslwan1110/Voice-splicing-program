% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   ���ܣ���Ƶ��ʽת���������ʽ��ת��           % 
% %   ���ߣ�Mr-Ma Technology(��ά)             %
% %   ʱ�䣺2020.04.09                           %
% %   ת����ע������                             %
% %   https://blog.csdn.net/qq_29225913/article/details/105445028
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% ֧�ֵ�ת����ʽ  %%%%%%%%%%%
% % |-----------------------------------------------------|
% % |            ��ʽ                     |  ���� |  ���  | 
% % |-----------------------------------------------------|
% % |WAVE (.wav)                          |  ֧�� |  ֧��  |
% % |OGG (.ogg)                           |  ֧�� |  ֧��  |
% % |FLAC (.flac)                         |  ֧�� |  ֧��  |
% % |AU (.au)                             |  ֧�� | ��֧�� |
% % |AIFF(.aiff/.aif)                     |  ֧�� | ��֧�� |
% % |AIFC (.aifc)                         |  ֧�� | ��֧�� |
% % |MP3 (.mp3) (ע����԰汾)             |  ֧�� | ��֧�� |
% % |MPEG-4 AAC(.m4a/.mp4)(ע����԰汾)   |  ֧�� |  ֧��  |
% % |-----------------------------------------------------|
% 
% % ������ϸ��Ϣ����ο� Matlab ���������ͨ�����������ʡ�λ�
% % audioread ��https://ww2.mathworks.cn/help/matlab/ref/audioread.html
% % audiowrite ��https://ww2.mathworks.cn/help/matlab/ref/audiowrite.html?ue&s_tid=gn_loc_dropp
% 
% %% ����ڴ���ı���  %%%%%%%%%%%
% clear; 
% close all; 
% clc;
% %%  ��������  %%%%%%%%%%%
% % ��Ƶ�ļ�����
% % source_FS = 16000; %ԭ�����ʣ���ֵ���ã����ļ���ȡ��
% dest_FS = 48000; % Ŀ�������
% bps = 16;        % λ�� -audiowrite ��������������� .wav/.flac �ļ���
% qlt = 80;       % Quality -audiowrite ��������������� .ogg �ļ���
% br  = 0;         % BitRate -audiowrite ��������������� .m4a/.mp4 �ļ��� 
% % Tfile = 0;       % ����Ƭ�ε���ʱ��s (0Ϊ���޸�) ��ʱû���������
% 
% % ��Ƶͨ��
% stereo = 'stereo';            % �����Ƶͨ�� 'mono'����������'stereo'��������
% audio_ch = 'left';          % ���������룬���������ʱ��ѡ����һͨ���� 'left'��'right'��'both'
% Left_channel_ratio = 0.5;   % ������ռ��
% 
% % �ļ�
% inputfile = ('Recording');   % �����ļ���
% outputfile = [inputfile,'.mp3'];    % ����ļ���
% AudioInfo = audioinfo(inputfile)    % ��ӡ������Ƶ�ļ�����Ϣ�����ӷֺţ�
% 
% %��������
% isplay = 1; % 0������ԭ��Ƶ 1������������Ƶ
% 
% % ��ͼ
% isdraw = 1;                  % 1���򿪻�ͼ�� 0���رջ�ͼ�� ���� FFT ��Ҫ����ʱ��
% t_axis = [-1 1];             % ʱ��ͼ��ʾ����(������)����ʵ������
% f_axis = [0 20000 0 0.002];  % Ƶ��ͼ��ʾ���򣬰�ʵ������
% 
% %%  ԭʼ��Ƶ�ļ���ȡ  %%%%%%%%%%%
% [y_input,source_FS]=audioread(inputfile);         % �� WAV �ļ�ת���ɱ���
% fs_str = ['(Fs:',num2str(source_FS/1000),'KHz)']; % fs ���ַ�������ͼ����ʾ�ļ� fs
% whos y_input                                      % ��ӡy����Ϣ
% if isplay == 0                                    % ���ѡ�񲥷�ԭ��Ƶ
%     if source_FS < 192000                         %     ����С�� 192K �����ʲ��ܲ���
%         sound(y_input,source_FS)                  %     ���ÿ��Բ��������ĺ��� sound()
%     end
% end
% 
%     
% %%  ��ӡ���벨��ͼ��  %%%%%%%%%%%
% if isdraw == 1
%     figure('Name','�����ļ�����','NumberTitle','off');
%     [N,n] = size(y_input);                  % �����ļ��ĳ�������������
%     if n==1                                 % ������
%         y1 = y_input(:,1);                  % ��ȡ��������
%         subplot(2,1,1);             
%         plot(y1);                           % ���������ͼ��
%         axis([0 N t_axis]);                 % ����ʵ����������ƿ���ʾ����
%         xlabel('������');ylabel('����');
%         title(['������ԭʼ���� ',fs_str]);
% 
%         subplot(2,1,2);
%         NFFT = source_FS * floor(N / source_FS);   % FFT ����ȡ fs �ı���
%         Y_fft = abs(fft(y1,NFFT));                 % ������ FFT ����
%         A_FFT1 = Y_fft * 2 / NFFT;                 % ת��Ϊ��ֵ
%         A_FFT1(1) = 0;                             % ����ֱ��������Ӱ��
%         f = (0:(NFFT/2)) * source_FS / NFFT;       % Ƶ�ʷ�Χ,�����꣨����ף�
%         plot(f ,A_FFT1(1:NFFT/2 +1));              % ��� FFT ͼ��
%         axis(f_axis);                              % ����ʵ����������ƿ���ʾ����
%         xlabel('Ƶ��/Hz');ylabel('���');
%         title('������ԭʼ����Ƶ��')
% 
%     else                                    % ������                          
%         y1 = y_input(:,1);                  % ��ȡ����������
%         y2 = y_input(:,2);                  % ��ȡ����������
% 
%         subplot(2,2,1);                     
%         plot(y1);                            % ���������ͼ��
%         axis([0 N t_axis]);                  % ����ʵ����������ƿ���ʾ����
%         xlabel('������');ylabel('����');
%         title(['������ԭʼ���� ',fs_str]);
% 
%         subplot(2,2,2);
%         plot(y2);                            % ���������ͼ��
%         axis([0 N t_axis]);                  % ����ʵ����������ƿ���ʾ����
%         xlabel('������');ylabel('����');
%         title(['������ԭʼ���� ',fs_str]);
% 
%         subplot(2,2,3); 
%         NFFT = source_FS * floor(N / source_FS);    % FFT ����ȡ fs �ı���
%         Y_fft1 = abs(fft(y1,NFFT));                 % ������ FFT ����
%         A_FFT1 = Y_fft1 * 2 / NFFT;                 % ת��Ϊ��ֵ
%         A_FFT1(1) = 0;                              % ����ֱ��������Ӱ��
%         f = (0:(NFFT/2)) * source_FS / NFFT;        % Ƶ�ʷ�Χ,�����꣨����ף�
% 
%         plot(f ,A_FFT1(1:NFFT/2 +1));               % ��� FFT ͼ��
%         axis(f_axis);                               % ����ʵ����������ƿ���ʾ����
%         xlabel('Ƶ��/Hz');ylabel('���');
%         title('������ԭʼ����Ƶ��')
% 
%         subplot(2,2,4);   
%         NFFT = source_FS * floor(N / source_FS);    % FFT ����ȡ fs �ı���
%         Y_fft2 = abs(fft(y2,NFFT));                 % ������ FFT ����
%         A_FFT2 = Y_fft2 * 2 / NFFT;                 % ת��Ϊ��ֵ
%         A_FFT2(1) = 0;                              % ����ֱ��������Ӱ��
%         f = (0:(NFFT/2)) * source_FS / NFFT;        % Ƶ�ʷ�Χ,�����꣨����ף�
% 
%         plot(f ,A_FFT2(1:NFFT/2 +1));               % ��� FFT ͼ��
%         axis(f_axis);                               % ����ʵ����������ƿ���ʾ����
%         xlabel('Ƶ��/Hz');ylabel('���');
%         title('������ԭʼ����Ƶ��')
%     end    
% end
% %% �����������  %%%%%%%%%%%
% if strcmp(stereo,'mono')                % ���ѡ��Ϊ������
%     if n == 1                           % ��������ļ��ǵ������ļ�
%         y_ratio = y1;                   % ��ô����ļ��������������ļ�������һ��
%     else                                % ����������������ļ�
%         if strcmp(audio_ch,'left')      %     ���ͨ��ѡ�� ������
%             y_ratio = y1;               %     ����������������Ϊ�������
%         elseif strcmp(audio_ch,'right') %     ���ͨ��ѡ�� ������
%             y_ratio = y2;               %     ����������������Ϊ�������
%         else                            %     ���ͨ��ѡ������ѡ�both��
%             y_ratio = y1*Left_channel_ratio + y2*(1-Left_channel_ratio); % �ͽ����������������������ݣ������趨�ı�������ϡ�
%         end
%     end
% else                                    % ������ѡ��Ϊ������
%     if n == 1                           % �����ļ��ǵ�����
%         y_ratio(:,1) = y1;              % ��ô���������������������
%         y_ratio(:,2) = y1;
%     else                                % �����ļ���������
%         y_ratio = y_input;              % ��ô����ļ��������������ļ�������һ��
%     end
% end
% 
% %% ���������Ƶ��ת��  %%%%%%%%%%%
% y_out = resample(y_ratio,dest_FS,source_FS);              % �Լ�����������ݣ������ز������趨�Ĳ�����
% audiowrite(outputfile,y_out,dest_FS,'BitsPerSample',bps); % ��������ݴ洢Ϊ��Ƶ�ļ���
% %audiowrite(outputfile,y_out,dest_FS,'Quality',qlt); % ��������ݴ洢Ϊ��Ƶ�ļ���
% 
% %% �����Ƶ�ļ���ȡ  %%%%%%%%%%%
% [y_output,dest_FS]=audioread(outputfile);            % ����һ�����ɵ��ļ���ȡ��������Ϊ����
% fs_str = ['(Fs:',num2str(dest_FS/1000),'KHz)'];      % �������ַ���
% 
% if isplay == 1                                       % ���ѡ�񲥷�������Ƶ
%     if dest_FS < 192000                              %     ����С�� 192K �����ʲ��ܲ���
%         sound(y_output,dest_FS)                      %     ���ÿ��Բ��������ĺ��� sound()
%     end
% end
% 
% %%  ��ӡ�������ͼ�� ������ӡ���벨��ͼ�����ƣ����ں��������Խϵͣ��Ͳ�������д�����ˣ�  %%%%%%%%%%%
% if isdraw == 1
%     figure('Name','����ļ�����','NumberTitle','off');    % �򿪵�2��figure
%     [N,n] = size(y_output);                    % �����ļ��ĳ������������� 
%     if n==1                                    % ������
%         y1 = y_output(:,1);                    % ��ȡ��������
% 
%         subplot(2,1,1);
%         plot(y1);                              % ���������ͼ��
%         axis([0 N t_axis]);                    % ����ʵ����������ƿ���ʾ����
%         xlabel('������');ylabel('����');
%         title(['���������ɲ��� ',fs_str]);
% 
%         subplot(2,1,2);
%         NFFT = dest_FS * floor(N / dest_FS);    % FFT ����ȡ fs �ı���
%         Y_fft = abs(fft(y1,NFFT));              % ������ FFT ����
%         A_FFT1 = Y_fft * 2 / NFFT;              % ת��Ϊ��ֵ
%         A_FFT1(1) = 0;                          % ����ֱ��������Ӱ��
%         f = (0:(NFFT/2)) * dest_FS / NFFT;      % Ƶ�ʷ�Χ,�����꣨����ף�
%         plot(f ,A_FFT1(1:NFFT/2 +1));           % ��� FFT ͼ��
%         axis(f_axis);                           % ����ʵ����������ƿ���ʾ����
%         xlabel('Ƶ��/Hz');ylabel('���');
%         title('���������ɲ���Ƶ��');
% 
%     else                                        % ������
%         y1 = y_output(:,1);                     % ��ȡ����������
%         y2 = y_output(:,2);                     % ��ȡ����������
% 
%         subplot(2,2,1);
%         plot(y1);                               % ���������ͼ��
%         axis([0 N t_axis]);                     % ����ʵ����������ƿ���ʾ����
%         xlabel('������');ylabel('����');
%         title(['���������ɲ��� ',fs_str]);
% 
%         subplot(2,2,2);
%         plot(y2);                               % ���������ͼ��
%         axis([0 N t_axis]);                     % ����ʵ����������ƿ���ʾ����
%         xlabel('������');ylabel('����');
%         title(['���������ɲ��� ',fs_str]);
% 
%         subplot(2,2,3);
%         NFFT = dest_FS * floor(N / dest_FS);    % FFT ����ȡ fs �ı���
%         Y_fft1 = abs(fft(y1,NFFT));             % ������ FFT ����
%         A_FFT1 = Y_fft1 * 2 / NFFT;             % ת��Ϊ��ֵ
%         A_FFT1(1) = 0;                          % ����ֱ��������Ӱ��
%         f = (0:(NFFT/2)) * dest_FS / NFFT;      % Ƶ�ʷ�Χ,�����꣨����ף�
% 
%         plot(f ,A_FFT1(1:NFFT/2 +1));           % ��� FFT ͼ��
%         axis(f_axis);                           % ����ʵ����������ƿ���ʾ����
%         xlabel('Ƶ��/Hz');ylabel('���');
%         title('���������ɲ���Ƶ��')
% 
%         subplot(2,2,4);                         
%         NFFT = dest_FS * floor(N / dest_FS);    % FFT ����ȡ fs �ı���
%         Y_fft2 = abs(fft(y2,NFFT));             % ������ FFT ����
%         A_FFT2 = Y_fft2 * 2 / NFFT;             % ת��Ϊ��ֵ
%         A_FFT2(1) = 0;                          % ����ֱ��������Ӱ��
%         f = (0:(NFFT/2)) * dest_FS / NFFT;      % Ƶ�ʷ�Χ,�����꣨����ף�
% 
%         plot(f ,A_FFT2(1:NFFT/2 +1));           % ��� FFT ͼ��
%         axis(f_axis);                           % ����ʵ����������ƿ���ʾ����
%         xlabel('Ƶ��/Hz');ylabel('���');
%         title('���������ɲ���Ƶ��')
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
