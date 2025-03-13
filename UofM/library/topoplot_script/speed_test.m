coord = rand(100,2);
val = rand(1,100);

% old = topoplot_general(val, coord,'shiftpreset', 3, 'smooth', 10);
new = topoplot_general_test(val, coord,'shiftpreset', 3, 'smooth', 10);