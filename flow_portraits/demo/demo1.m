% Nathaniel Linden
% Bing Brunton Lab - University of Washington
% June 29th, 202

% This script provides a tutorial for creating a simple FLOW portrait for the data stored in Demo_Data.mat
% The data is provided by Professor William Moody (Biology, University of Washington)

clear all; close all; clc
addpath('../flow_portraits')

% load data
load('../../Data/DemoData.mat')

% Example Usage of creating a flow portrait
% Here we specify the filename and the inegration length
% We also enable the white_background function for the final image to set the area around the cortex to white

filename        = './example1';
integration_len = 40;
flowPortrait(dat.dfof, integration_len, 'white_background', true, 'mask', dat.mask, 'filename', filename);

