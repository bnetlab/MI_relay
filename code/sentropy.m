% calculate shannon entropy of a discrete pdf
function H = sentropy(p)
H = sum(-(p(p>0).*(log2(p(p>0)))));