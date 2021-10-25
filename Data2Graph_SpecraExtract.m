%%muplti_cau_disperision_Ref_G150_C680_50X
%%本版本针对2021.10.23_角分辨数据_G150
%%Part1:修改了光谱的采集范围 129.52-224.5
%%Part2:加入了针对反射谱设置的归一化计算
%%Part3:使用了 50 倍镜头
%%Part4:添加了提取Abs谱信息

clc;
clear;
close all; 
file=dir('F:\matlab_project\2021.10.15_ARR_PVSK_FIB_Pattern\Data-txt\*.txt');
n_max = length(file);
y = []; y(1) = 0; y(2) = 173;
AA= 96; BB= -0.0618;   % y = a + b * x 边缘作直线拟合得到倾斜角度关系

Grating = 600 ; Center = 530;    % 使用光谱仪采谱的设定参数，光栅刻线（Grating）和中心波长（Center）
xx_res = 5000;  yy_res = 1000;   % 画图插值使用的分辨率，xx_res对应波长轴，yy_res对应角度轴
NA = 0.8; n_air = 1;             % 镜头数值孔径NA，空气折射率n_air

  fr_bg = importdata('90013-BG-G600-C530-1s-5t_1.txt');  % 导入背景文件数据
% fr_Lamp = importdata('9001-Lamp-50X-g150-C620-500ms-10t_1.txt');   % 导入灯谱文件数据 9011-Lamp-50X-g150-C620-500ms-10t-P274_1
  fr_Lamp = importdata('99923-Ag-G600-C530-1s-5t_1.txt');   % 导入灯谱文件数据 9011-Lamp-50X-g150-C620-500ms-10t-P274_1
 
for n =1:n_max
  %strposition = strfind(file(n).name, '_g150_c570'); %仅导入反射光谱的数据
  strposition = strfind(file(n).name, '23012-S23-G600-C530-1s-5t-P056_1');
  if ~isempty (strposition)
  fr_dep = importdata(file(n).name);  % 导入主文件数据
  
  %% 光谱倾斜修正
    fr_dep_Line2_2 = fr_dep(:,2) - (fr_dep(:,1) .* BB)-AA;  %对应波长修正数据
    %fr_dep_fixed = [fr_dep_Line1 fr_dep_Line2_2 fr_dep_Line3];
    fr_dep = [fr_dep(:,1) fr_dep_Line2_2 fr_dep(:,3)];  
    % 此时的矩阵已经是矫正倾斜后的数据
  
  %% 减去背景，除以白光灯谱
      % fr_dep_2 = fr_dep(:,3).*10 - fr_bg(:,3);
      fr_dep_2 = fr_dep(:,3) - fr_bg(:,3);
      fr_Lamp_2 = fr_Lamp(:,3) - fr_bg(:,3);
      fr_dep_2 = fr_dep_2 ./ fr_Lamp_2;
      fr_dep = [fr_dep(:,1:2) fr_dep_2]; % 重新组成3矩阵形式  % 此时的矩阵已经是减去仪器噪声背景（并除以灯谱白光）后的数据

  %% 对角度进行转换修正
      %fr_dep_fixed(:,3)= log10(fr_dep(:,3));%设置强度显示为对数,不用时设置为注释
      [k]=find(fr_dep(:,2)>= round( 0.5*(y(1)+y(2))- 0.5*abs((y(2)-y(1))))&fr_dep(:,2)<= round(0.5*(y(1)+y(2))+ 0.5*abs((y(2)-y(1)))));
      fr_dep_region = fr_dep(k,1:3); % 筛选符合 y(1)-y(2)区间的数据
      arg = atand(linspace(-tand(asind(NA/n_air)),tand(asind(NA/n_air)),length(fr_dep_region))); %构建角度坐标轴 ？？
      arg_max = asind(NA/n_air);   %角度坐标轴的最大值
      arg_min = -arg_max;
      % % % % 构建新的y轴坐标
      fr_dep_region(:,2) = fr_dep_region(:,2) - (0.5*(y(2)-y(1))); % 数据矩阵依照y轴中心线对齐
      %L_max = (0.5*(y(2)-y(1)));
      L0 = (0.5*(y(2)-y(1)));                                   % 计算采集角对应变长，为后续角度换算作准备
      fr_dep_region(:,2) = asind(fr_dep_region(:,2)./ L0);       % 对数据矩阵进行角度左边的换算操作
      %fr_dep_region(:,2) = arg;  % % ？？？ 这玩意儿插入替换就废了
      xx_region = linspace(max(fr_dep_region(:,1)),min(fr_dep_region(:,1)),xx_res); %定义x轴坐标
      %yy_region = linspace(min(fr_dep_region(:,2)),max(fr_dep_region(:,2)),yy_res); %定义y轴坐标
      yy_region = linspace(arg_min,arg_max,yy_res); %定义y轴坐标   %%%% 改变了选取y坐标的范围
      yy_region = yy_region';
      for jj = 1:yy_res
         xx_region(jj,:) = xx_region(1,:);
      end
      for jj = 1:xx_res
        yy_region(:,jj) = yy_region(:,1);
      end
      zz_region = griddata(fr_dep_region(:,1),fr_dep_region(:,2),fr_dep_region(:,3),xx_region,yy_region);
      energy_region = zeros(size(xx_region));
      energy_region = energy_region+1240;
      xx_region = energy_region./xx_region;

  %% 画图
      figure
      pcolor(yy_region,xx_region,zz_region);shading interp%伪彩色图
      set(gcf,'Colormap',bone)
      xlim([-arg_max arg_max])
     
      if (1)
      if Grating==150 & Center==640    % Grating150_Center640nm
          ylim([1.52 2.56])
      elseif Grating==150 & Center==840    % Grating150_Center840nm
          ylim([1.22 1.86])    
      elseif Grating==150 & Center==570    % Grating150_Center570nm
          ylim([1.665 2.57])
      elseif Grating==600 & Center==530    % Grating600_Center530nm
          ylim([2.167 2.52])
      elseif Grating==600 & Center==540    % Grating600_Center530nm
          ylim([2.13 2.49])
      end  
      % ylim([1.52 2.56]) %[Grating150_Center640nm] (464.561, 814.497) (2.669eV,1.522eV)
      % ylim([2.167 2.52]) %[Grating600_Center530nm] (487.345, 572.289) (2.544eV,2.167eV)
      % ylim([1.22 1.86]) %[Grating150_Center840nm] (665, 1013) (1.86eV,1.22eV)
      end
      
      xlabel('Angle (degree)','FontName','Times newman','FontSize',10);
      ylabel('Energy (eV)','FontName','Times newman','FontSize',10);%设置xy轴标签内容和字体
      set(gca, 'Fontname', 'Times newman', 'Fontsize', 10);%设置xy轴的字体类型和字号大小的
      set(gca,'Position',[.18 .17 .6 .73]);%这句是设置xy轴在图片中占的比例，可能需要自己微调。
       Mesh_size = xx_res * yy_res ;
       zz_region_1 = reshape(zz_region,1,Mesh_size);

      Zmax= max(zz_region_1);
      Zmin= min(zz_region_1);
      
      if Zmax >1 
          Cmax = 1; %设置color bar 显示范围
      else Cmax = Zmax*0.95; % Cmax=Zmax-0.03; %设置color bar 显示范围
      end
      if Zmin <0 
          Cmin=0; %设置color bar 显示范围
      else Cmin= Zmin*1.03; % Cmin= Zmin+0.02; %设置color bar 显示范围
      end
       set(gca,'CLim',[Cmin Cmax]); %设置color bar 显示范围 
       
    if (1)
      if     Grating==150 & Center==640    % Grating150_Center640nm
          set(gca,'Ytick',[1.5:0.2:2.6])
      elseif Grating==150 & Center==840    % Grating150_Center840nm
          set(gca,'Ytick',[1.2:0.1:1.8])     
      elseif Grating==150 & Center==570    % Grating150_Center570nm
          set(gca,'Ytick',[1.4:0.2:2.6])
      elseif Grating==600 & Center==530    % Grating600_Center530nm
          set(gca,'Ytick',[2.2:0.1:2.5])
      elseif Grating==600 & Center==540    % Grating600_Center530nm
          set(gca,'Ytick',[2.2:0.1:2.5])
      end 
    end
      
      % set(gca,'CLim',[0.005 0.1]); %设置color bar 显示范围  
      set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);%这句是将线宽改为2
      % set(gca,'Ytick',[1.86:0.3:2.9]);%设置xy轴的刻度显示步长 %[Grating150_Center550nm]
      % set(gca,'Ytick',[2.2:0.1:2.5]);
      % set(gca,'Ytick',[1.2 :0.1:1.8]); % g150 c840

      set(gca,'Xtick',[-50:25:50]);
      set(gcf,'Position',[100 100 300 300]);%这句是设置绘图的大小，不需要到word里再调整大小。我给的参数，图的大小是7cm
      colorbar('position',[0.82 0.17 0.03 0.73]); %设置彩色条的位置和大小

       %print([num2str(file(n).name),'-MCD','.tif'] ,'-dtiffn','-r300');
      
  %%提取Abs谱信息
         Abs = [0, 0.5, 1, 2, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52];
             [~,Length_Abs]=size(Abs);
             fr_y_0 = yy_region(:,1);
             Spectra=[];
        if(1)
                 for i = 1:1:Length_Abs
                     Abs_i = Abs(1,i);
                     [~,Abs_i_1]=min(abs(fr_y_0(:) - Abs_i)); 
                     Spec_i = zz_region(Abs_i_1,:); 
                     Spec_i = [Abs_i, fr_y_0(Abs_i_1,1), Spec_i];
                     Spectra = [Spectra;Spec_i];  
                 end
        end

             fr_x_0 = xx_region(1,:);
             Wavelength = [0,0,fr_x_0];
             Spectra = [Wavelength;Spectra];
             Spectra_1 = Spectra';
           %   save(file(n).name,'Spectra_1','-ascii');
        

  end

end
