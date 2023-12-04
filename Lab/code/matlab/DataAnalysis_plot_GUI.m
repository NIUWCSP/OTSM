load('ScattorData.mat');
load('BERData.mat');
AllBERDataAverage = sum(AllBERData')./length(AllBERData(1,:));

% Create a new figure
fig = figure('Position', [100 100 800 400]);

% 創建複選框
num_lines = length(AllRxDataSymbEqAverage(1,:))/18;
checkbox_width = 100;
checkbox_height = 20;
checkbox_spacing = 30;
checkbox_top = 20;
checkbox_left = 20;
for i = 1:num_lines
    checkbox_label = sprintf('Line %d', mod(-i,8));
    checkbox_position = [checkbox_left, checkbox_top + (i-1)*checkbox_spacing, checkbox_width, checkbox_height];
    uicontrol('Style', 'checkbox', 'String', checkbox_label, 'Value', 1, 'Position', checkbox_position, 'Callback', {@checkbox_callback, i});
end
            
% Store the checkbox and button handles in the figure's UserData
set(fig, 'UserData', struct('cb', cb));

% Create subplots for the scatter plot and equalized scatter plot
subplot(121);
plot(0:length(AllBERDataAverage(1,:))-1,AllBERDataAverage);axis([-0.1,length(AllBERDataAverage(1,:))-1+0.1,0,1]);title('BPSK vs. OFDM')

subplot(122);legend_text = cell(1, length(AllRxDataSymbEqAverage(1,:))/18);

for i = 1:length(AllRxDataSymbEqAverage(1,:))/18
    plot(reshape(AllRxDataSymbEqAverage(1:108,18*(i-1)+1:18*i),1,108*18).*exp(-1i*pi/4),'.');
    axis([-1.5,1.5,-1.5,1.5]);
    title('Scatter plot after equalization');
    axis square;
    hold on;
    legend_text{i}=[' Distance = ',num2str(i-1),' cm '];
end
% % 預設將1-7cm關閉(勾選單問題未解決)
% scotter_h = get(gca, 'Children');
% for i = 2:length(AllRxDataSymbEqAverage(1,:))/18-1
%     set(scotter_h(i), 'Visible', 'off');
% end

legend(legend_text);

%% 關閉顯示圖形
function checkbox_callback(hObject, eventdata, line_index)
    scotter_h = get(gca, 'Children');
    if get(hObject, 'Value')
        set(scotter_h(line_index), 'Visible', 'on');
    else
        set(scotter_h(line_index), 'Visible', 'off');
    end
end