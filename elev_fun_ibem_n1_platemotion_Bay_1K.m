function elev_pred = elev_fun_ibem_n1_platemotion_Bay_1K(m, G_tect, G_bay, d_outlets, Ginv_elev, channel_elements, compare_elements)

% number of observation points
Nx = length(G_tect);
 
% unknown parameters
 
tectslip = m(1:2);
n = 1;
 
% F = nan(Nx,2);
% Ngeo = length(geo_map(:,1)); 
% logKcon     = m(3:(2+Ngeo));
logKcon     = m(3);

F = G_tect ./ (10.^logKcon);

% for i = 1:Ngeo
%     F(geo_map(i,2):geo_map(i,3),:) = G_tect(geo_map(i,2):geo_map(i,3),:) ./ (10.^logKcon(geo_map(i,4)));
% end
 
ks = (F*tectslip).^(1./n);

d = [ks;d_outlets];

uplift_constrain=(G_bay*tectslip)*1e3; % in [mm/yr]
% SAfar_constrain=(G_SAfar*tectslip)*1e4;

elev_pred = Ginv_elev * d;  % in [m]

elev_pred = elev_pred(channel_elements);
elev_pred = elev_pred(compare_elements);

% elev_pred=[elev_pred;uplift_constrain;SAfar_constrain];
elev_pred=[elev_pred;uplift_constrain];

