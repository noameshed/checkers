function new_block = otsu_local_thresh(block)
% This function uses Otsu's algorithm to threshold the image regionally

T = graythresh(block);
new_block = block<255*T;


