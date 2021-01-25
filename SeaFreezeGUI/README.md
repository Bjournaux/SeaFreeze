# SeaFreeze GUI (Coming by early february 2021)
The SeaFreeze GUI allows thermodynamic and seismic data to be effortlessly viewed and extracted for ices Ih, II, III, V, and VI and liquid water.
## Installation 
Select the corresponding version for your operating system. Using the app installer, download the standalone app. Once loaded, the app can be opened and run. All necessary files are included within the package. 
## Using the Interface 
The interface includes a series of panels and tab subsets for exploration of thermodynamic and seismic properties. 
![Screenshot 2020-12-01 171323](https://user-images.githubusercontent.com/70338087/100821704-c89ac880-3405-11eb-94eb-d47613ee0f94.jpg)
### Step 1: Select a Material and Units
To visualize and export data, select a material and choose the desired units.
### Step 2: Set Parameters and Plot EOS 
Type the desired pressure (MPa) and temperature ranges (K) into min/max field boxes. The chosen values should reflect the valid range included next to the material name. Input the number of points in the nP/nT boxes. Next, click the 'Plot EOS' button; the data is now ready for export. 
*Note: Data must be plotted and ranges must be set correctly for exportation.* 
### Step 3:Export and Explore 
All generated ranged data can be exported in the lower left-hand corner of the screen. Check the desired options, select a format, and click 'Export'. 
The SeaFreeze GUI also includes a 'Single Point Calculator' and 'Water Phase Diagram' to examine and extract data. 

### Reading the Exported Files 
#### Excel files:
##### Metadata
A metadata sheet is printed in addition to other data sheets. This sheet includes the version of SeaFreeze, date printed, 
a second row containing min, max, and number of points for temperature, followed a third row that houses min, max, and the number of points for pressure. The number of points is used to create a linearly spaced vector for temperature and pressure ranges.  

![Metadataexcel](https://user-images.githubusercontent.com/70338087/93945188-48f8cb00-fceb-11ea-8198-e41b602c6764.png)


##### Thermodynamic and Seismic Data 
All generated data is printed on separate sheets labeled by the sheet name (property (units)). The data matrix presents in the following way: 

The number of rows correspond to the linearly spaced pressure points and the number of columns refer to the evenly spread temperature coordinates; the first row and first column showcases the minimum pressure and temperature respectively.
i.e if 200 MPa is the minimum pressure, 400 MPa is max pressure, and the number of points is 101, the first row represents a pressure of 200 MPa, the second row showcases 220 MPa etc.... 

![2020-09-22 (8)](https://user-images.githubusercontent.com/70338087/93942773-ee10a500-fce5-11ea-9a94-671950cfc200.png)

#### Txt and Json Files
##### Metadata
Metadata reads the same as 'Excel Files', but prints above the outputted data and under the name/units of the property. 
##### Thermodynamic and Seismic Data 
The data formats the same as the 'Excel Files'. However, if there are too many columns, the data prints into 'groups'. These sets showcase one complete row. 

![txtdata2](https://user-images.githubusercontent.com/70338087/93945547-ff5cb000-fceb-11ea-9e94-3f84631c15a1.png)

#### Water Phase Diagram Excel and Txt Files 
##### Metadata 
The date and version of SeaFreeze prints at the beginning of both file types.
##### Pressure and Temperature Coordinate Data 
Two columns appear below the metadata. The first column showcases temperature points in Kelvin. Its paired pressure (MPa) coordinate is adjacent in column two.
If triple points are included in export, the layout in Excel column orients as above, except it includes an extra column on the right to indicate materials.
The triple points for txt files have sets of 3 lines: materials in line 1, temperature coordinates in line 2, and pressure coordinates in line 3. i.e. (material;T;P). 

![txtwpd](https://user-images.githubusercontent.com/70338087/93951068-3043e180-fcfa-11ea-8891-92be53f272a3.png)

##### Changes in Enthalpy and Entropy
###### MetaData
Both excel and text metadata includes date and version along with first and last paired temperature/pressure points. The number of points used to generate the temperature and pressure domain follows the pairs.
*The identity of the materials can be found in either the file name (.txt) or the sheet name (.xlsx).*
###### Reading the Change in Enthalpy and Entropy
The column on the left belongs to the change in enthalpy of a given phase equilibrium. The second column is the change of entropy on the chosen phase line. The first listed value for the change in enthalpy and entropy corresponds to the first metadata temperature and pressure coordinate. Delta H and S values end on the last listed temperature and pressure point in the metadata, where points in between relate to linearly spaced temperature and pressure values.

##### Excel

![Excel HS](https://user-images.githubusercontent.com/70338087/100824257-b1120e80-340a-11eb-9f8b-cd8604a3581e.png)

##### Text

![TEXT HS](https://user-images.githubusercontent.com/70338087/100824247-ae171e00-340a-11eb-8018-3fb37e4db9c0.png)

## Examples

### Visually Displaying and Exporting Data 
To export data (Specific Heat, Density, and Internal Energy) for ice Ih: 
Material 'Ice Ih' and units are selected. The pressure and temperature min/max are changed to reflect the valid range of 'Ice Ih'. 
The 'Plot EOS' button is clicked:
The desired data is checked on the lower left-hand corner and the export button is clicked.
The data is saved to a user-specified directory.

![Examples](https://user-images.githubusercontent.com/70338087/100825515-123ae180-340d-11eb-8a2d-ad0964fdd1a4.png)

### Single Point Calculator 
Values for ice Ih at 0.1 MPa and 273.1 K in J/mol are calculated in the following manner; 
'Ih' and units are selected from the materials list under 'Graphical Options'.
Calculate is selected.
The properties and potentials for the material are now displayed and can be copied and pasted. 

![Single Point](https://user-images.githubusercontent.com/70338087/100824240-ab1c2d80-340a-11eb-8fd4-091cfcee920d.png)

### Water Phase Diagram 
To generate the full Sea Freeze water phase diagram, check the appropriate box. 
Click on the 'Plot' button to visualize the data. 
The line coordinates are saved by choosing a file type and exporting. 

![WPD](https://user-images.githubusercontent.com/70338087/100824270-b8391c80-340a-11eb-8f8d-384fea445a7a.png)

## Remarks 
Range of validity
SeaFreeze stability prediction is currently considered valid down to 130K, which correspond to the ice VI - ice XV transition. The ice Ih - II transition is potentially valid down to 73.4 K (ice Ih - ice XI transition).
## Authors 
**Erica Clinton**- *University of Washington, Earth and Space Sciences Department, Seattle, USA* 

**Baptiste Journaux** - *University of Washington, Earth and Space Sciences Department, Seattle, USA*

