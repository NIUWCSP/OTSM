function [H,P]= Gen_delay_sequency_channel_matrix(N,M,G,Wn)
P=zeros(N*M,N*M);                  
for j=1:N
    for i=1:M
        E=zeros(M,N);
        E(i,j)=1;
        P((j-1)*M+1:j*M,(i-1)*N+1:i*N)=E; %equation (35) in [R3]
    end
end
H=kron(eye(M),Wn)*(P'*G*P)*kron(eye(M),Wn); %equation (34) in [R1]
end