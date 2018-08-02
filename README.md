# Image Processing for Board Game Piece Recognition

    Noam Eshed
    May 3, 2018
    ECSE 4540 Intro to Digital Image Processing

## Summary
This project used the MatLab Image Processing Toolbox to detect, recognize, and locate checker pieces on boards.
The goal is to detect all of the game pieces on the board, and return their locations (i.e. in the  board below, pieces are located at A1, A3, A5, etc.)

![alt text](https://github.com/Gnome42/checkers/blob/master/Report%20Images/labeled_checkerboard.jpg)

## Process Overview

### Board Transformation
This program begins by transforming the board image into a square. The user is prompted to select the four corners of the board, and a geometric transformation is applied accordingly.

### Canny Edges
The program first finds all of the edges on the chess board, then continues to find the pieces within them.
I found that using the image's red channel provided the most accurate edges in later processing.
Canny edge detection is applied to the red channel to find the edges in the board, then morphological opening is applied to clean up noise.

![alt text](https://github.com/Gnome42/checkers/blob/master/Report%20Images/canny_edges.JPG)

### Hough Edges
The Hough transform is used next to extract the strongest line edges in the image. 
In the next step, these edges are extended to the edges of the image, and improbable edges (i.e. those which are not near vertical or horizontal) are removed.
The Hough transform does not usually find all of the board edges, so mathematical interpolation is done to find the remaining edges.

##### Raw Hough Edges
![alt text](https://github.com/Gnome42/checkers/blob/master/Report%20Images/raw_hough.JPG)

##### Interpolated Edges
![alt text](https://github.com/Gnome42/checkers/blob/master/Report%20Images/all_edges.jpg)

### Thresholding for Game Pieces
I used 6 thresholding methods to find the checker pieces. Each method has its own strengths and weaknesses, and combined, they were able to detect most or all pieces on the board.
The methods are:
1. Original (no thresholding)
2. Red Channel (grayscale)
3. Global Otsu thresholding
4. Regional Otsu Thresholding
5. Simple Thresholding ( im_red_channel < 55 )
6. Canny Edges

MatLab function imfindcircles was applied to the result of each of these thresholding methods. The function looked for circles of radii between 40 and 60. Since the image transformation at the beginning of the entire process creates a 1000x1000 pixel image, this is relatively consistent.

![alt text](https://github.com/Gnome42/checkers/blob/master/Report%20Images/thresholds.jpg)

### Game Piece Processing
Since the six methods often produced multiple circles for each piece, the script then filters through all of the saved game piece locations and averages those with similar centroid coordinates.
While this is imperfect, it reduces most excess circles on the board. The result shows the board game with edges and pieces marked, and prints out the board location of each piece:

In this case, the resulting image is below and the listed checker locations are:
    A2, A4, A6, A8, B1, B3, B5, B7, C2, C4, C4, C6, C6, C6, C8, F1, F3, F5, F7, G2, G4, G4, G6, G6, G6, H1, H3, H5, H7
    
![alt text](https://github.com/Gnome42/checkers/blob/master/Report%20Images/final.jpg)


For more detail about the image processing algorithms and my process, see my final report [here](https://github.com/Gnome42/checkers/blob/master/Final%20Report.pdf).