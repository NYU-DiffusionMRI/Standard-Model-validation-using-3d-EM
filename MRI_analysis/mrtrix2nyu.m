function dki_nyu = mrtrix2nyu(tensor_mrtrix,kurtosis_mrtrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NYU (dki_fit)
% % Tensor
% D_ind =
%      1     1
%      1     2
%      1     3
%      2     2
%      2     3
%      3     3
% % Kurtosis
% W_ind =
%1.       1     1     1     1
%2.       1     1     1     2
%3.       1     1     1     3
%4.       1     1     2     2
%5.       1     1     2     3
%6.       1     1     3     3
%7.       1     2     2     2
%8.       1     2     2     3
%9.       1     2     3     3
%10.      1     3     3     3
%11.      2     2     2     2
%12.      2     2     2     3
%13.      2     2     3     3
%14       2     3     3     3
%15       3     3     3     3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % MRtrix     
% % Tensor
% volumes 0-5: D11, D22, D33, D12, D13, D23
% % Kurtosis
% volumes 0-2: W1111., W2222., W3333 
% volumes 3-8: W1112., W1113., W1222., W1333., W2223., W2333. 
% volumes 9-11: W1122., W1133., W2233. 
% volumes 12-14: W1123., W1223., W1233.
% % MRtrix (1st index 0 -> 1)   
% % Tensor
% volumes 1-6: D11, D22, D33, D12, D13, D23
% % Kurtosis
% volumes 1-3: W1111., W2222., W3333 
% volumes 4-9: W1112., W1113., W1222., W1333., W2223., W2333. 
% volumes 10-12: W1122., W1133., W2233. 
% volumes 13-15: W1123., W1223., W1233.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
st = size(tensor_mrtrix);
sk = size(kurtosis_mrtrix);
if(~isequal(st(1:3),sk(1:3)))
error('tensor and kurtosis should have the first three spatial sizes equal')
end
% dti: nyu->mrtrix = [1 4 6 2 3 5]
% dki: nyu->mrtrix = [1 11 15 2 3 7 10 12 14 4 6 13 5 8 9]

% dti: mrtrix->nyu = [1 4 5 2 6 3]
% dki: mrtrix->nyu = [1 4 5 10 13 11 6 14 15 7 2 8 12 9 3]
t2      = tensor_mrtrix(:,:,:,[1 4 5 2 6 3]); 
k2      = kurtosis_mrtrix(:,:,:,[1 4 5 10 13 11 6 14 15 7 2 8 12 9 3]);
dki_nyu = cat(4,t2,k2);
end % 