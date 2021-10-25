%%muplti_cau_disperision_Ref_G150_C680_50X
%%���汾���2021.10.23_�Ƿֱ�����_G150
%%Part1:�޸��˹��׵Ĳɼ���Χ 129.52-224.5
%%Part2:��������Է��������õĹ�һ������
%%Part3:ʹ���� 50 ����ͷ
%%Part4:�������ȡAbs����Ϣ

clc;
clear;
close all; 
file=dir('F:\matlab_project\2021.10.15_ARR_PVSK_FIB_Pattern\Data-txt\*.txt');
n_max = length(file);
y = []; y(1) = 0; y(2) = 173;
AA= 96; BB= -0.0618;   % y = a + b * x ��Ե��ֱ����ϵõ���б�Ƕȹ�ϵ

Grating = 600 ; Center = 530;    % ʹ�ù����ǲ��׵��趨��������դ���ߣ�Grating�������Ĳ�����Center��
xx_res = 5000;  yy_res = 1000;   % ��ͼ��ֵʹ�õķֱ��ʣ�xx_res��Ӧ�����ᣬyy_res��Ӧ�Ƕ���
NA = 0.8; n_air = 1;             % ��ͷ��ֵ�׾�NA������������n_air

  fr_bg = importdata('90013-BG-G600-C530-1s-5t_1.txt');  % ���뱳���ļ�����
% fr_Lamp = importdata('9001-Lamp-50X-g150-C620-500ms-10t_1.txt');   % ��������ļ����� 9011-Lamp-50X-g150-C620-500ms-10t-P274_1
  fr_Lamp = importdata('99923-Ag-G600-C530-1s-5t_1.txt');   % ��������ļ����� 9011-Lamp-50X-g150-C620-500ms-10t-P274_1
 
for n =1:n_max
  %strposition = strfind(file(n).name, '_g150_c570'); %�����뷴����׵�����
  strposition = strfind(file(n).name, '23012-S23-G600-C530-1s-5t-P056_1');
  if ~isempty (strposition)
  fr_dep = importdata(file(n).name);  % �������ļ�����
  
  %% ������б����
    fr_dep_Line2_2 = fr_dep(:,2) - (fr_dep(:,1) .* BB)-AA;  %��Ӧ������������
    %fr_dep_fixed = [fr_dep_Line1 fr_dep_Line2_2 fr_dep_Line3];
    fr_dep = [fr_dep(:,1) fr_dep_Line2_2 fr_dep(:,3)];  
    % ��ʱ�ľ����Ѿ��ǽ�����б�������
  
  %% ��ȥ���������԰׹����
      % fr_dep_2 = fr_dep(:,3).*10 - fr_bg(:,3);
      fr_dep_2 = fr_dep(:,3) - fr_bg(:,3);
      fr_Lamp_2 = fr_Lamp(:,3) - fr_bg(:,3);
      fr_dep_2 = fr_dep_2 ./ fr_Lamp_2;
      fr_dep = [fr_dep(:,1:2) fr_dep_2]; % �������3������ʽ  % ��ʱ�ľ����Ѿ��Ǽ�ȥ�������������������Ե��װ׹⣩�������

  %% �ԽǶȽ���ת������
      %fr_dep_fixed(:,3)= log10(fr_dep(:,3));%����ǿ����ʾΪ����,����ʱ����Ϊע��
      [k]=find(fr_dep(:,2)>= round( 0.5*(y(1)+y(2))- 0.5*abs((y(2)-y(1))))&fr_dep(:,2)<= round(0.5*(y(1)+y(2))+ 0.5*abs((y(2)-y(1)))));
      fr_dep_region = fr_dep(k,1:3); % ɸѡ���� y(1)-y(2)���������
      arg = atand(linspace(-tand(asind(NA/n_air)),tand(asind(NA/n_air)),length(fr_dep_region))); %�����Ƕ������� ����
      arg_max = asind(NA/n_air);   %�Ƕ�����������ֵ
      arg_min = -arg_max;
      % % % % �����µ�y������
      fr_dep_region(:,2) = fr_dep_region(:,2) - (0.5*(y(2)-y(1))); % ���ݾ�������y�������߶���
      %L_max = (0.5*(y(2)-y(1)));
      L0 = (0.5*(y(2)-y(1)));                                   % ����ɼ��Ƕ�Ӧ�䳤��Ϊ�����ǶȻ�����׼��
      fr_dep_region(:,2) = asind(fr_dep_region(:,2)./ L0);       % �����ݾ�����нǶ���ߵĻ������
      %fr_dep_region(:,2) = arg;  % % ������ ������������滻�ͷ���
      xx_region = linspace(max(fr_dep_region(:,1)),min(fr_dep_region(:,1)),xx_res); %����x������
      %yy_region = linspace(min(fr_dep_region(:,2)),max(fr_dep_region(:,2)),yy_res); %����y������
      yy_region = linspace(arg_min,arg_max,yy_res); %����y������   %%%% �ı���ѡȡy����ķ�Χ
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

  %% ��ͼ
      figure
      pcolor(yy_region,xx_region,zz_region);shading interp%α��ɫͼ
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
      ylabel('Energy (eV)','FontName','Times newman','FontSize',10);%����xy���ǩ���ݺ�����
      set(gca, 'Fontname', 'Times newman', 'Fontsize', 10);%����xy����������ͺ��ֺŴ�С��
      set(gca,'Position',[.18 .17 .6 .73]);%���������xy����ͼƬ��ռ�ı�����������Ҫ�Լ�΢����
       Mesh_size = xx_res * yy_res ;
       zz_region_1 = reshape(zz_region,1,Mesh_size);

      Zmax= max(zz_region_1);
      Zmin= min(zz_region_1);
      
      if Zmax >1 
          Cmax = 1; %����color bar ��ʾ��Χ
      else Cmax = Zmax*0.95; % Cmax=Zmax-0.03; %����color bar ��ʾ��Χ
      end
      if Zmin <0 
          Cmin=0; %����color bar ��ʾ��Χ
      else Cmin= Zmin*1.03; % Cmin= Zmin+0.02; %����color bar ��ʾ��Χ
      end
       set(gca,'CLim',[Cmin Cmax]); %����color bar ��ʾ��Χ 
       
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
      
      % set(gca,'CLim',[0.005 0.1]); %����color bar ��ʾ��Χ  
      set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);%����ǽ��߿��Ϊ2
      % set(gca,'Ytick',[1.86:0.3:2.9]);%����xy��Ŀ̶���ʾ���� %[Grating150_Center550nm]
      % set(gca,'Ytick',[2.2:0.1:2.5]);
      % set(gca,'Ytick',[1.2 :0.1:1.8]); % g150 c840

      set(gca,'Xtick',[-50:25:50]);
      set(gcf,'Position',[100 100 300 300]);%��������û�ͼ�Ĵ�С������Ҫ��word���ٵ�����С���Ҹ��Ĳ�����ͼ�Ĵ�С��7cm
      colorbar('position',[0.82 0.17 0.03 0.73]); %���ò�ɫ����λ�úʹ�С

       %print([num2str(file(n).name),'-MCD','.tif'] ,'-dtiffn','-r300');
      
  %%��ȡAbs����Ϣ
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
