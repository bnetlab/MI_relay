% calculate differential entropy of a continous pdf
function H = diff_entropy(p,bin_size)
H = bin_size*sum(-(p(p>0).*(log2(p(p>0)))));