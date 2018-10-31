# Matlab Algorithm for Neuron Detection from Images
Matlab algorithm and user interface for neuron network extraction from images

## User Interface
![Graphical interface screenshot](https://raw.githubusercontent.com/valeriu-balaban/matlab-neuron-detection-2d/master/user-interface.png)

1. Open tif stack - Opens the select file dialog window. After loading the image, the program searches for .mat file with the same name as the image to load the previously saved location of neurons.
2. Tif frame selection slider. Moreover, the left and right arrow keys move the slider one position to the left or to the right. 
3. The image containing neurons for the detection. Clicking on an empty region will triggers the semi-automatic detection algorithm that searches to find the radius for the neuron with the center at the click location. A yellow square appears to indicate the borders of the neuron. A subsequent click on the square would select it and change its color from yellow to magenta. Pressing 'escape' removes the selection. To change the position of neuron, select it and the click on an empty region to move its center to the click location. To delete a neuron, press the D key while a neuron is selected. 
4. The detection algorithm. The algorithm searches for the ring of white pixels that often surrounds neurons and based on the ratio of how much it surrounds the region, that region is selected as a neuron or not.
  * Neuron Radius a list of radii in the Matlab format (start-value:step-increment:end-value).
  * Edge Width specifies the width of the white ring that surround the neurons.
  * Theta Threshold sets how much at least the ring should surround the region for it to be selected as a neuron. E.g., a value of 1 means the ring should completely surround the regions, and 0.5 the ring should be a semicircle.
5. Manual adjustment of the radius. The slider is enable only when a neuron is selected. Changing the position of the scrollbar would change the neuron radius.
6. Location of all neurons and their respective radius.

All changes to the neurons are automatically saved to the disk.

__ATTENTION__ The undo functionality is not implemented.  To avoid unintentional selection or repositioning press often the escape key.
