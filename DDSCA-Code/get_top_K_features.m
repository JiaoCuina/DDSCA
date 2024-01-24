function idx_selected = get_top_K_features(w, top_k_selected)

if nargin < 2 
	top_k_selected = 10;
end

absw = abs(w);
[w_sorted, w_sorted_idx] = sort(absw, 'descend');
idx_selected = w_sorted_idx(1 : max(top_k_selected));