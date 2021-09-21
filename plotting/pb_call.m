%% Subfunctions for buttons
% function pb_call2(varargin)
% %     disp('NOT Flipped');
%     uiresume;
%     return
% end

function pb_call(hObject,~)
    hObject.UserData = struct('moveon',1);
%     disp('Flipped');
    uiresume;
    return
end