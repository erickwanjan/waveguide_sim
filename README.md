Subwavelength Grating (SWG) Waveguide Simulator

Part of the research project "Analyzing Subwavelength Gratings for Low-Loss Waveguide Design" done by Eric Jan. Code and simulator designed to predict the effect of using sbwavelength grating waveguides on loss. Intended to be used for future low loss waveguide design for fabrication by bulk CMOS process.

All code is contained in the utils.py and utils_client.py files. For descriptions of functionality, consult the code documentation of each file.

DEPENDENCIES
• python (compatible with both python2.7 and python3)
• matplotlib, particularly matplotlib.pyplot and matplotlib.patches (for plotting)
• numpy

HOW TO RUN:

	• METHOD 1
		simply run the utils_client.py file interactively
		eg. run:
		-------------------------
		python -i utils_client.py
		-------------------------

	• METHOD 2
	Step 1
		import utils_client.py into your python file
		eg. run:
		-------
		python3
		-------
	Step 2
		in the python environment, run the following
		---------------------------
		from utils_client import *
		import utils_client
		---------------------------
