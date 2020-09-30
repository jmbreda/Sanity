function [mi,ari,fm,jaccard] = get_mi_ari_fm_jaccard(my_class,my_cluster)
	% get normalized mutual information
	mi = nmi(my_class,my_cluster);

    % Get Confision Matrix " cm(1,1) = TP, cm(1,2) = FN, cm(2,1) = FP, cm(2,2) = TN
	cm = zeros(2,2);
	N_cell = length(K);
	for i = 1:(N_cell-1)
		for j = (i+1):N_cell
			Truth = G(i)==G(j);
			Seen  = K(i)==K(j);
			% Confusion Matrix
			% cm(1,1) = TP
			% cm(1,2) = FN
			% cm(2,1) = FP
			% cm(2,2) = TN
			cm(-Truth+2,-Seen+2) = cm(-Truth+2,-Seen+2)+1;
		end
	end

	% Compute ARI, FM and Jaccard
    ari = 2*(cm(1,1)*cm(2,2) - cm(1,2)*cm(2,1))/( (cm(1,1)+cm(1,2))*(cm(1,2)+cm(2,2)) + (cm(1,1)+cm(2,1))*(cm(2,1)+cm(2,2)) );
    fm  = sqrt( (cm(1,1)/(cm(1,1)+cm(1,2))) * (cm(1,1)/(cm(1,  1)+  cm(2,1))) );
    jaccard = cm(1,1)/(cm(1,1)+cm(1,2)+cm(2,1));


