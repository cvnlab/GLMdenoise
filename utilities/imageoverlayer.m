function imageoverlayer(underlay, overlay, thresh)
    
    % TODO - make some guess as to if volumes need to be skipped. 
    skip = 1;

    if size(underlay,3) >1

        ax1 = axes;
        imagesc(makeimagestack(underlay(:,:,1:skip:end)));
        colormap(ax1, 'gray')
        axis('image');
        hold on
        % lets mask some things, to see stuff. 
        alphahelper = overlay(:,:,1:skip:end)>thresh;
        alphahelper = alphahelper + (overlay(:,:,1:skip:end)<(-1*thresh));
        ax2 = axes;
        imagesc(ax2, makeimagestack((overlay(:,:,1:skip:end))), 'AlphaData', makeimagestack(alphahelper))
        axis('image')
        colormap(ax2, 'jet')
        ax2.Visible = 'off';
        
    else
        % its just a slice, get crazy
        ax1 = axes;
        imagesc(underlay);
        colormap(ax1, 'gray')
        axis('image')
        hold on
        % lets mask some things, to see stuff. 
        alphahelper = overlay>thresh;
        alphahelper = alphahelper + (overlay<(-1*thresh));
        ax2 = axes;
        imagesc(ax2, overlay, 'AlphaData', alphahelper)
        axis('image')
        colormap(ax2, 'jet')
        ax2.Visible = 'off';
    end
