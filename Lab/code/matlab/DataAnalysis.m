load('ScattorData.mat');
load('BERData.mat');
AllBERDataAverage = sum(AllBERData')./length(AllBERData(1,:));
subplot(121);plot(0:length(AllBERDataAverage(1,:))-1,AllBERDataAverage);axis([-0.1,length(AllBERDataAverage(1,:))-1+0.1,0,1]);title('錯誤率對距離圖');
subplot(122);
for i = 1:length(AllRxDataSymbEqAverage(1,:))/18-1:length(AllRxDataSymbEqAverage(1,:))/18
    plot(reshape(AllRxDataSymbEqAverage(1:108,18*(i-1)+1:18*i),1,108*18).*exp(-1i*pi/4),'.');axis([-1.5,1.5,-1.5,1.5]);title('scattor after equalization'); axis square;
    hold on;
    legend_text{i}=[' 距離 = ',num2str(i-1),' cm '];
    
end

legend(legend_text{1},legend_text{8});

% for i = 1:length(AllRxDataSymbEqAverage(1,:))/18
%     plot(reshape(AllRxDataSymbEqAverage(1:108,18*(i-1)+1:18*i),1,108*18).*exp(-1i*pi/4),'.');axis([-1.5,1.5,-1.5,1.5]);title('scattor after equalization'); axis square;
%     hold on;
%     legend_text{i}=[' 距離 = ',num2str(i-1),' cm '];
%     
% end
% legend(legend_text);

set(gcf, 'Position', get(0, 'Screensize'));