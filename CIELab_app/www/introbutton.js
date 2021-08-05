document.addEventListener('DOMContentLoaded', function() {
	var selector = document.getElementById("matrix");
		selector.setAttribute("data-step", "1");
		selector.setAttribute("data-intro", "First step: Choose the type of data you will upload (or test by using the examples). For the purpose of this tutorial, please choose any of these tracks. 'Single-cell' option exprects a matrix counts with cells as rows and genes as columns. 'High-dimensional' refers to the rest of the matrices that are high dimensional but they will not fo through 'Single-cell RNAseq analysis'. 'Distance matrix' requests a matrix with equal number of rows and columns and '3D data' can be any dataset with 3 numeric columns and optional a 4th column with names that users want to be translated into colors.");
		selector.setAttribute("data-position", "bottom");
		
		selector = document.getElementById("browse.well");
		selector.setAttribute("data-step", "2");
		selector.setAttribute("data-intro", "For the purposes of this tutorial, please choose the given example.");
		selector.setAttribute("data-position", "bottom");
		
		selector = document.getElementById("contents");
		selector.setAttribute("data-step", "3");
		selector.setAttribute("data-intro", "Check the first 10 rows of your dataset.");
		selector.setAttribute("data-position", "bottom");
		
		selector = document.getElementById("btnanalysis");
		selector.setAttribute("data-step", "4");
		selector.setAttribute("data-intro", "'Analysis' button is now activated and needs to be pressed to start your analysis.");
		selector.setAttribute("data-position", "bottom");
		
		selector = document.getElementById("table_box");
		selector.setAttribute("data-step", "5");
		selector.setAttribute("data-intro", "The optimization algorithm gives more than one solutions. In the 3D plot you can see the most optimal solution, while from this table you can choose another one. The main difference will be the colors but you can clearly understand how the point cloud is spread in the tab named '2D projections.'");
		selector.setAttribute("data-position", "bottom");
		
		selector = document.getElementById("sliders");
		selector.setAttribute("data-step", "6");
		selector.setAttribute("data-intro", "Teasing the sliders you can change the direction that the point cloud is spread and as a result the colors. The optimization algorithm will find the 'best' solution based on the 'new' axes. Lstar refers to the brightness, Astar changes the range between blue-yellow and Bstar between green-red. Finally, you can increase the scaling factor in order to have more intense colors.");
		selector.setAttribute("data-position", "right");
		
		selector = document.getElementById("weightButton");
		selector.setAttribute("data-step", "7");
		selector.setAttribute("data-intro", "After changing the sliders of the axes, press this button to re-optimize based on your new choices.");
		selector.setAttribute("data-position", "bottom");
		
		selector = document.getElementById("scalingButton");
		selector.setAttribute("data-step", "8");
		selector.setAttribute("data-intro", "After changing the scaling factor, press this button to re-optimize based on your new choices.");
		selector.setAttribute("data-position", "bottom");
		
		selector = document.getElementById("satellite_div");
		selector.setAttribute("data-step", "9");
		selector.setAttribute("data-intro", "Here you can check the 2D colored view of your dataset and how the point cloud is spread across the axes (especially after re-optimizing with new values in sliders).");
		selector.setAttribute("data-position", "bottom");
		
	  selector = document.getElementById("satellite");
		selector.setAttribute("data-step", "10");
		selector.setAttribute("data-intro", "The plots are interactive. By clicking and dragging your cursor over the points you can select them and check the colored legend. Deselect with double-click on an empty space.");
		selector.setAttribute("data-position", "bottom");

	  selector = document.getElementById("legend");
		selector.setAttribute("data-step", "11");
		selector.setAttribute("data-intro", "Color legend. If you select points on the cloud, you will be able to see only the selected genes with their colors.");
		selector.setAttribute("data-position", "bottom");
		
	  selector = document.getElementById("download_div");
		selector.setAttribute("data-step", "12");
		selector.setAttribute("data-intro", "Download your tsv file here! A preview is shown below with the table.");
		selector.setAttribute("data-position", "bottom");
		
	  selector = document.getElementById("openModal");
		selector.setAttribute("data-step", "13");
		selector.setAttribute("data-intro", "For further information, questions or bug reports, find contact information here. Enjoy!");
		selector.setAttribute("data-position", "bottom");

	
	document.getElementById("introButton").onclick = function(){
		introJs().start();
	};
}, false);
