function ave_hours = fun_hours_means(pol_hours,stat_dist)
%Calculate average hours for each quintile q1,..,q5

[n_a,n_z] = size(pol_hours);

if ~isequal(size(stat_dist),size(pol_hours))
    error('pol_hours and stat_dist have different shapes!')
end

pol_hours = reshape(pol_hours,[n_a*n_z,1]);
stat_dist = reshape(stat_dist,[n_a*n_z,1]);

[~,sort_ind] = sort(pol_hours);
pol_hours = pol_hours(sort_ind);
stat_dist = stat_dist(sort_ind);
cum_stat_dist = cumsum(stat_dist);

[~,p20_cut] = min(abs(cum_stat_dist-0.2));
[~,p40_cut] = min(abs(cum_stat_dist-0.4));
[~,p60_cut] = min(abs(cum_stat_dist-0.6));
[~,p80_cut] = min(abs(cum_stat_dist-0.8));

% h_p20 = pol_hours(p20_cut);
% h_p40 = pol_hours(p40_cut);
% h_p60 = pol_hours(p60_cut);
% h_p80 = pol_hours(p80_cut);

ave_hours = zeros(5,1);
mass      = zeros(5,1);

q1_ind = 1:p20_cut;
q2_ind = p20_cut+1:p40_cut;
q3_ind = p40_cut+1:p60_cut;
q4_ind = p60_cut+1:p80_cut;
q5_ind = p80_cut+1:numel(stat_dist);


ave_hours(1) = sum(pol_hours(q1_ind).*stat_dist(q1_ind));
ave_hours(2) = sum(pol_hours(q2_ind).*stat_dist(q2_ind));
ave_hours(3) = sum(pol_hours(q3_ind).*stat_dist(q3_ind));
ave_hours(4) = sum(pol_hours(q4_ind).*stat_dist(q4_ind));
ave_hours(5) = sum(pol_hours(q5_ind).*stat_dist(q5_ind));

mass(1) = sum(stat_dist(q1_ind));
mass(2) = sum(stat_dist(q2_ind));
mass(3) = sum(stat_dist(q3_ind));
mass(4) = sum(stat_dist(q4_ind));
mass(5) = sum(stat_dist(q5_ind));

% n_x = n_a*n_z;
% for ix = 1:n_x
%     if pol_hours(ix)< h_p20
%         mass(1) = mass(1)+stat_dist(ix);
%         ave_hours(1) = ave_hours(1)+pol_hours(ix)*stat_dist(ix);
%     elseif pol_hours(ix)>=h_p20 && pol_hours(ix)<h_p40
%         mass(2) = mass(2)+stat_dist(ix);
%         ave_hours(2) = ave_hours(2)+pol_hours(ix)*stat_dist(ix);
%     elseif pol_hours(ix)>=h_p40 && pol_hours(ix)<h_p60
%         mass(3) = mass(3)+stat_dist(ix);
%         ave_hours(3) = ave_hours(3)+pol_hours(ix)*stat_dist(ix);
%     elseif pol_hours(ix)>=h_p60 && pol_hours(ix)<h_p80
%         mass(4) = mass(4)+stat_dist(ix);
%         ave_hours(4) = ave_hours(4)+pol_hours(ix)*stat_dist(ix);
%     else 
%         mass(5) = mass(5)+stat_dist(ix);
%         ave_hours(5) = ave_hours(5)+pol_hours(ix)*stat_dist(ix);    
%     end %end if
% end %end ix

ave_hours = ave_hours./mass;

end %end function