function [bp_sum, bp_sum_cf] = kirchhoff_forward_3d(R01, d_01, R0, R1, d_R0, d_R1, range_profiles, index_mat, phase_mat)
c=3.0e8;
%% BP³ÉÏñ
echo = range_profiles;
tr_pair_number = size(range_profiles, 2);
[a,b] = size(range_profiles);
range_profiles_1 = zeros(a,b);
range_profiles_2 = zeros(a,b);
%{
for i = 2:a
    range_profiles(i,:) = echo(a-i+2,:);
end
%}

for m = 1:tr_pair_number
    range_profiles_1(1:a-1,m)=diff(range_profiles(:,m));
    range_profiles_2(1:a-2, m)= diff(range_profiles(:,m),2);
    
end

channel_images = zeros(size(index_mat));
for i = 1:tr_pair_number
    temp = range_profiles(:, i);
    temp_1 = range_profiles_1(:,i);
    temp_2 = range_profiles_2(:,i);
    %
    channel_images(:, :, :, i) = phase_mat(:,:,:,i).*d_01(:,:,:,i).*R01(:,:,:,i).*...
                                 ((temp_2(index_mat(:,:,:,i))/(c^2))+((1./R0(:,:,:,i)+1./R1(:,:,:,i))/c.*temp_1(index_mat(:,:,:,i)))...
                                 +(temp(index_mat(:,:,:,i))./R01(:,:,:,i)));%匹配
       %}     
    
%{
    channel_images(:, :, :, i) = 4*phase_mat(:,:,:,i).*d_01(:,:,:,i).*...
            (R01(:,:,:,i)/(c^2).*temp_2(index_mat(:,:,:,i))+(R0(:,:,:,i)+...
            R1(:,:,:,i))/c.*temp_1(index_mat(:,:,:,i))+temp(index_mat(:,:,:,i)));%Zhuge
%}
    %{
channel_images(:, :, :, i) = 4*phase_mat(:,:,:,i).*d_01(:,:,:,i).*...
            (R01(:,:,:,i)/(c^2).*temp_2(index_mat(:,:,:,i))+(R0(:,:,:,i)+...
            R1(:,:,:,i))/c.*d_R0(:,:,:,i).*temp(index_mat(:,:,:,i))+temp(index_mat(:,:,:,i)));%Liu Xin
    %}
end

% Ô­Ê¼BP
bp_sum = abs(sum(channel_images, 4));

% CF¼ÓÈ¨BP
cf = abs(bp_sum).^2 ./ (tr_pair_number*sum(abs(channel_images).^2, 4));
bp_sum_cf = bp_sum.*cf;