
%%%% testing the current code %%%%
clc;
clear all;
close all;




BW_MyMethod2 = im2bw(imread('Binary/dna-3/MyMethod2.png'));
BW_GT = im2bw(imread('Binary/dna-3/GT.png'));
BW_CP = im2bw(imread('Binary/dna-3/CP.png'));
[RI_CP,JI_CP] = RandJaccardIndex(BW_GT,BW_CP)
[RI_MM,JI_MM] = RandJaccardIndex(BW_GT,BW_MyMethod2)







