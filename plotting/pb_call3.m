%% Subfunctions for buttons
function pb_call3(hObject,~)
    disp('EXIT');
    hObject.Parent.UserData = 'EXIT';
    uiresume;
    return
end