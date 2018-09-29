# Matlab Algorithm for Neuron Detection from Images
Matlab algorithm and user interface for neuron network extraction from images

## User Interface
![Graphical interface screenshot](https://raw.githubusercontent.com/valeriu-balaban/matlab-neuron-detection-2d/master/user-interface.png)

## Algorithm description

1. Open tif stack - Shows the open file dialog window allowing selection of files with .tif and .tiff extensions. After the program reads all the images from the tif stack it tries to load the previous save neurons location from a .mat file with the same name as the image it if exists.
2. Tif frame selection slider. Also, the left and right arrow keys trigger the transition to the previous and respectively next frame.
3. The image of neurons. A click on an empty region will trigger the detection algorithm for a neuron with the center at the click location. After completion, a yellow square would appear. A subsequent click on the square would select it, and its color would change from yellow to magenta. Next, any clicks on empty regions would move the neuron to that position, or change the selection if the click was inside another neuron. To delete a neuron, press the D key while a neuron is selected. Pressing 'escape' removes the selection.
4. The detection algorithm. The algorithm searches for the ring of white pixels that commonly surrounds the neurons and based on the ratio of how much it surrounds the regions is selected as a neuron or not.
  * Neuron Radius defines the radii used by detection algorithm in the Matlab list format (start-value:step-increment:end-value).
  * Edge Width specifies the width of the white rings that surround the neurons.
  * Theta Threshold sets the lowest ratio of the edge ring surrounding the neurons. E.g., a value of 1 means the rings completely surrounds the neurons, and for 0.5 the ring forms a semicircle around the neuron.
5. Manual adjustment of the radius. Enabled only if a neuron is selected and changing the position of the scrollbar would change the neuron radius.
6. Location of all neurons and their respective radius.

All changes to the neurons are automatically saved to the disk.

__ATTENTION__ The undo functionality is not implemented.  To avoid unintentional selection or repositioning press often the escape key.
