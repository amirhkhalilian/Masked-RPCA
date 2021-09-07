function Visual_view3d(img,rec)

switch nargin
    case 1
        [m,n,z] = size(img);
        mini = min(img(:));
        maxi = max(img(:));
        slice = 1;

        fig1 = figure;
        im = imagesc(img(:,:,slice)); axis image;
        colormap('gray'); colorbar; caxis([mini maxi]);

        drawnow;

        txt = uicontrol(fig1,...
                'Style','text',...
                'String','slice 1',...
                'Units', 'Normalized',...
                'Position', [0.2 0.05 0.6 0.05]);
        sld = uicontrol(fig1,...
                'Style', 'slider',...
                'Min',1,'Max',z,'Value',1,...
                'Units', 'Normalized',...
                'Position', [0.2 0.0 0.6 0.05],...
                'Callback', {@fun1, img, txt, im },...
                'SliderStep', [1/(z-1) 1]);
    case 2
        if ~isreal(rec)
            warning('showing the real part of the recovered image');
        end
        [m,n,z] = size(img);
        slice = 1;

        fig1 = figure;
        ax(1)=subplot(1,2,1); im1 = imagesc(img(:,:,slice)); colormap gray; axis image; title('image 1');caxis([0,1]);
        ax(2)=subplot(1,2,2); im2 = imagesc(real(rec(:,:,slice))); colormap gray; axis image; title('image 2');caxis([0,1]);
        linkaxes(ax,'xy');

        drawnow;

        txt = uicontrol(fig1,...
                'Style','text',...
                'String','slice 1',...
                'Units', 'Normalized',...
                'Position', [0.2 0.05 0.6 0.05]);
        sld = uicontrol(fig1,...
                'Style', 'slider',...
                'Min',1,'Max',z,'Value',1,...
                'Units', 'Normalized',...
                'Position', [0.2 0.0 0.6 0.05],...
                'Callback', {@fun2, img, rec, txt, im1, im2 },...
                'SliderStep', [1/(z-1) 1]);
    otherwise
        error('not enough input arguments');
end

end

function fun1(hObject, eventdata, img, txt, im)
slice = floor(get(hObject,'Value'));
im.CData = img(:,:,slice);
txt.String = ['slice ',num2str(slice)];
end

function fun2(hObject, eventdata, img, rec, txt, im1, im2)
slice = floor(get(hObject,'Value'));
im1.CData = img(:,:,slice);
im2.CData = rec(:,:,slice);
txt.String = ['slice ',num2str(slice)];
end