function [img] = mergePatch(p, y, x, Y, X)
  img = zeros(Y, X);
  coeff = zeros(Y, X);
  p_idx = 1;
  for xx=1:x
      for yy=1:y
          pp = col2im(p(p_idx,:), [y x], [Y X], 'sliding');
          img(yy:yy+Y-y,xx:xx+X-x) = img(yy:yy+Y-y,xx:xx+X-x)+pp;
          coeff(yy:yy+Y-y,xx:xx+X-x) = coeff(yy:yy+Y-y,xx:xx+X-x)+1;
          p_idx = p_idx+1;
      end
  end
  img = img ./ coeff;
  return;
end
