classdef FakeExtractFuture < handle
    %FAKEEXTRACTFUTURE Stands in for extractTrial's future in processResult.
    properties
        Predecessors
    end
    methods
        function obj = FakeExtractFuture()
            obj.Predecessors = FakeExtractFuture.empty;
        end
        function [H, S, LS, F0, SNR] = fetchOutputs(~)
            [H, S, LS, F0, SNR] = deal([]);
        end
    end
end
