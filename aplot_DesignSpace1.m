'''
Created by Junu Kim (j-kim@pse.t.u-tokyo.ac.jp) at the University of Tokyo
Modified on Jan. 25th, 2024

Publication title: Kinetic study and model-based design space determination for a drug substance flow synthesis using amination reaction via nucleophilic aromatic substitution
'''

for xxx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    filename1 = sprintf('1_Design_space_total_ID_%d.csv', xxx)
    %filename1='1_PMI_DP_total_f_1_ID_0_temp_' + xxx + '.csv';%プロットしたいデータのファイル指定
    a=csvread(filename1);%読込
    filename2='1_axis.csv';%プロットする軸データのファイル指定
    b=csvread(filename2);%読込
    
    c= [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]%b(find(b(:,1)),1); %x軸の定義
    d= [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0] %b(find(b(:,2)),2); %y軸の定義

    %a(a == 1) = 3
    %a(a == 0) = 1
    %a(a == 3) = 0
    
    contourf(c,d,a,1,'LineColor','k', 'LineWidth', 1.5)%等高線作成

    map=[179 179 179;
    255	255	255];
    
    map = map / 256.0
    
    colormap(map)
    %set(gca,'clim','FontSize',16)%軸の数値のフォントサイズ
    %set(gca,'clim',[0, 1],'FontSize',16)%軸の数値のフォントサイズ
    set(gca,'clim',[min(a,[],'all') , Inf ],'FontSize',16)%軸の数値のフォントサイズ
    %colorbar
    pbaspect([1 1 1])

    save_filename = sprintf('1_Design_space_total_ID_%d.tif', xxx)
    saveas(gcf, save_filename)
end

for xxx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    filename1 = sprintf('1_Design_space_total_PMI_RprodID_%d.csv', xxx)
    %filename1='1_PMI_DP_total_f_1_ID_0_temp_' + xxx + '.csv';%プロットしたいデータのファイル指定
    a=csvread(filename1);%読込
    filename2='1_axis.csv';%プロットする軸データのファイル指定
    b=csvread(filename2);%読込
    
    c= [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]%b(find(b(:,1)),1); %x軸の定義
    d= [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0] %b(find(b(:,2)),2); %y軸の定義

    %a(a == 1) = 3
    %a(a == 0) = 1
    %a(a == 3) = 0
    
    contourf(c,d,a,1,'LineColor','k', 'LineWidth', 1.5)%等高線作成

    map=[179 179 179;
    255	255	255];
    
    map = map / 256.0
    
    colormap(map)
    %set(gca,'clim','FontSize',16)%軸の数値のフォントサイズ
    %set(gca,'clim',[0, 1],'FontSize',16)%軸の数値のフォントサイズ
    set(gca,'clim',[min(a,[],'all') , Inf ],'FontSize',16)%軸の数値のフォントサイズ
    %colorbar
    pbaspect([1 1 1])

    save_filename = sprintf('1_Design_space_total_PMI_Rprod_ID_%d.tif', xxx)
    saveas(gcf, save_filename)
end
