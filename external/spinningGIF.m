function spinningGIF(fname)
% spinningGIF(fname): makes a spinning GIF of the current plot and saves it
% Usage: make your 3D plot (using plot3(...) or scatter3(...) etc.) and
% then call SpinningGIF with the file name that you want
% 
% All credit to Adina Stoica at
% https://uk.mathworks.com/matlabcentral/fileexchange/46050-spinning-gif

axis off
%     view(0,10)
center = get(gca, 'CameraTarget');
pos = get(gca, 'CameraPosition');
radius = norm(center(1:2) - pos(1:2));
angles = 0:0.01*pi:2*pi;
for ii=1:(length(angles)-1)
   angle = angles(ii);
   set(gca, 'CameraPosition', [center(1) + radius * cos(angle),...
                               center(2) + radius * sin(angle),...
                               pos(3)]);
   drawnow;
   frame = getframe(1);
   im = frame2im(frame);
   [imind,cm] = rgb2ind(im,256);
   if ii == 1
       imwrite(imind,cm,fname,'gif', 'Loopcount',inf);
   else
       imwrite(imind,cm,fname,'gif','WriteMode','append','DelayTime', 0.02);
   end
end


end