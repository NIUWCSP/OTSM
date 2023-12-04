% 打開第一個.fig文件
fig1 = openfig('BPSK.fig');
fig2 = openfig('OFDM.fig');

% 獲取第一個圖形的子圖
axes1 = findobj(fig1, 'Type', 'Axes');
axes2 = findobj(fig2, 'Type', 'Axes');

figure;
subplot(1,2,1);title('錯誤率對距離圖');
copyobj(allchild(axes1(2)), gca);
hold on;
copyobj(allchild(axes2(2)), gca);set(copyobj(allchild(axes2(2)), gca),'Color',[1 0 0]);legend('','BPSK','OFDM');
subplot(1,2,2);axis([-1.5,1.5,-1.5,1.5]);title('scatter after equalization'); axis square;
legend('距離 = 0 cm','距離 = 7 cm');
%copyobj(allchild(axes1(1)), gca);
hold on;
copyobj(allchild(axes2(1)), gca);

close(fig1);
close(fig2);

set(gcf, 'Position', get(0, 'Screensize'));