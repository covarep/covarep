function delta_mat = get_delta_mat(mat)

% Function to calculate deltas from a feature matrix.

delta_mat=zeros(size(mat));
N=size(mat,2);

for n=1:N
    delta_mat(2:end,n)=diff(mat(:,n));
end