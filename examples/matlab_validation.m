% Run this manually in MATLAB with the ORIGINAL AbsWinding.m on your path.
% Fixtures are clean, finite, with no repeated consecutive points.
fixtures = {'straight',[0 0;1 0;2 0]; 'right_angle',[0 0;1 0;1 1]; 'u_turn',[0 0;1 0;0 0]; 'square_open',[0 0;1 0;1 1;0 1;0 0]};
for k=1:size(fixtures,1)
    xy=fixtures{k,2};
    data=[(0:size(xy,1)-1)',xy];
    fprintf('%s,%.17g\n',fixtures{k,1},AbsWinding(data));
end
% Expected: straight 0; right_angle .25; u_turn .5; square_open .75.
% square_open has 3 interior turns; the original does NOT close tangent endpoints.
