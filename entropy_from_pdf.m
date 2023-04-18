function entropy = entropy_from_pdf(pdf)

entropy = nansum(arrayfun(@(x)(-x.*log2(x)), pdf),'all');