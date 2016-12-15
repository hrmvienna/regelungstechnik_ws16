function pd = polydiff(p)
    pd = []
    len = length(p)
    for k = 1:1:(len-1)
        pd = [pd, p(k) * (len - k)]
   end
end