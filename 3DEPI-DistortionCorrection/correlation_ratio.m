function eta = correlation_ratio(X, Y)
X = X(:); % make sure we've got column vectors, simplifies things below a bit
Y = Y(:);
L = max(X);
mYx = zeros(1, L+1); % we'll write mean per class here
nx = zeros(1, L+1);  % we'll write number of samples per class here
for i = unique(X).'
   Yn = Y(X == i);
   if numel(Yn)>1
      mYx(i+1) = mean(Yn);
      nx(i+1) = numel(Yn);
   end
end
mY = mean(Y);        % mean across all samples
eta = sqrt(sum(nx .* (mYx - mY).^2) / sum((Y-mY).^2));