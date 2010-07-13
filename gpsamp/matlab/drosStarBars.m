function drosStarBars(h, pvals),

% Move the texts slightly because Matlab misaligns them
xfudge = -.05;
for k=1:size(pvals, 2),
  x0 = get(get(h(k), 'Children'), 'XData');
  x = mean(x0(2:3, :)) + xfudge;
  y = get(h(k), 'YData');
  
  for l=1:length(x),
    if y(l) < 100 && y(l) > 0,
      if pvals(l, k) < .001,
	ht = text(x(l), y(l)+1, '***', 'Rotation', 90, 'VerticalAlignment', 'cap', 'HorizontalAlignment', 'left');
      elseif pvals(l, k) < .01,
	ht = text(x(l), y(l)+1, '**', 'Rotation', 90, 'VerticalAlignment', 'cap', 'HorizontalAlignment', 'left');
      elseif pvals(l, k) < .05,
	ht = text(x(l), y(l)+1, '*', 'Rotation', 90, 'VerticalAlignment', 'cap', 'HorizontalAlignment', 'left');
      end
    end
  end
end
