%% Subfunctions for buttons
function pb_call2(hObject,~)
    disp('');
    %At the end of the callback function:
    hObject.Parent.UserData = 'REMOVE';
    uiresume;
    return
end