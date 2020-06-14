# Culture Data Viewer
View culture data with pandas and plotly®


## What does it do
With Culture Data Viewer (CDV) you can work with culture data coming directly from the reactor\* and merge your offline measurements. CDV calculates automatically relevant parameters (e.g., yields, rates), and lets you explore beautiful interactive charts made with plotly.
<br>
![Culture Data Viewer in action](CDV2.PNG?raw=true)

### Input
\* Culture Data Viewer is currently developed to work with data from DASGIP® reactors.

CDV can load data coming from offline measurements - these data need to follow the structure defined in [Culture_Data_Viewer_example.ipynb](Culture_Data_Viewer_example.ipynb).


### Output
CDV lets you explore interactively the imported culture data in a jupyter notebook with plotly.


### Something more
CDV's calc module computes some of the most common parameters from your culture data: find the average yield of a compound, the instant production rate, the total amount of compound produced.


## Getting started


### Prerequisites
In order to use CDV you will need python with the pandas, numpy, plotly and plotly-express libraries.
CDV was tested to work with:
```
python == 3.8.2 
pandas == 1.0.3
numpy == 1.18.4
plotly == 4.8.0
plotly-express == 0.4.1
```

### Usage

CDV can be imported into a jupyter notebook. Here is an example of how to use it: [Culture_Data_Viewer_example.ipynb](Culture_Data_Viewer_example.ipynb)




## Status: work in progress!

Current limitations:
Culture Data Viewer was initially developed to work with data coming from DASGIP® reactors. For this reason only these specific files can be loaded.

If you find any bug or if you have any interesting ideas, feel free to contribute to the project! :)


 



## License

This project is licensed under the GNU LGPLv3 - see the [LICENSE](LICENSE) file for details.





<br><br><br>

### Notice of Non-Affiliation and Disclaimer

This software is not affiliated, associated, authorized, endorsed by, or in any way officially connected with Eppendorf, or any of its subsidiaries or its affiliates. The official Eppendorf website can be found at www.eppendorf.com. The name “DASGIP” as well as related names, marks, emblems and images are registered trademarks of Eppendorf. 
