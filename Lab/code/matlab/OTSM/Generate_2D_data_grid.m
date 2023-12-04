function X = Generate_2D_data_grid(N,M,x_data,data_grid)
x_vec=zeros(N*M,1);
data_array=reshape(data_grid,1,N*M);
[~,data_pos]=find(data_array>0);
x_vec(data_pos)=x_data;
X=reshape(x_vec,M,N);
end