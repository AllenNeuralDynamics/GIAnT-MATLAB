function tf = isInteractiveSession()
%ISINTERACTIVESESSION True if a person can respond to dialogs and figure windows.
%
% False when MATLAB was started with -batch (e.g. scheduled or pipeline runs),
% cannot show figure windows (-nodisplay, as on headless servers such as Code
% Ocean), or runs without Java. Pipeline steps use this to avoid opening a GUI
% that would wait forever for input.

tf = ~batchStartupOptionUsed && usejava('jvm') && feature('ShowFigureWindows');
end
