# CMS HH limit extraction and combination package

The purpose of this package is to have centralized CMS HH code base for the limit extraction and combination between the different channels for both model independent and model dependent cases.

This package is based on the tools, code and interfaces provided by [CombineHarvester](https://github.com/cms-analysis/CombineHarvester) package and [CMS Higgs Combination toolkit](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit).

## How to install

The installation steps follows guidelines from CombineHarvester (http://cms-analysis.github.io/CombineHarvester/) and CMS combine (https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit/) documentations.

The recommended way to install HHStatAnalysis framework is by using install.sh script located in https://github.com/cms-hh/HHStatAnalysis/blob/master/install.sh
```shell
curl -s https://raw.githubusercontent.com/cms-hh/HHStatAnalysis/master/install.sh | bash -s [MODE] [N_JOBS] [CMSSW_RELEASE]
```
where MODE is the installation mode: *full* for the full installation, and *plotting* to install only the plotting-related code (default *full*); N_JOBS is the number of jobs to run simultaneous during the compilation (default 8), CMSSW_RELEASE is the CMSSW release (default CMSSW_7_4_7).

Example:
```shell
curl -s https://raw.githubusercontent.com/cms-hh/HHStatAnalysis/master/install.sh | bash -s
```


## How to run

For details about how to run tools provided by CombineHarvester (e.g. combineTool.py, plotLimits.py), please, see its documentation: http://cms-analysis.github.io/CombineHarvester/.


### How to run limits for a single channel

- Create a configuration file (see the next section).
- Run all limit extraction steps up to plotting using run_hh_limits.py script:
```shell
run_hh_limits.py --cfg CFG_FILE --model-desc DESC --output OUTPUT_PATH [--parallel N_JOBS] shape_file [shape_file] ...
```
where CFG_FILE is the configuration file; DESC is the name of a model descriptor entry in the config; OUTPUT_PATH is the path where the limit files should be stored; N_JOSB is the number of parallel jobs; shape_file - file (or files) with the input shapes.

Examples:
```shell
# bbtautau resonant model independent limits
run_hh_limits.py --cfg HHStatAnalysis/Resources/LimitSetups/ttbb_resonant.cfg --model-desc ttbb_res --output output/limits shapes/LLR_shapes.root
# bbtautau resonant limits for hMSSM model
run_hh_limits.py --cfg HHStatAnalysis/Resources/LimitSetups/ttbb_resonant.cfg --model-desc ttbb_res_hMSSM --output output/limits_hMSSM shapes/LLR_shapes.root
```

## Configuration description

### Overview

Used text-based configuration format originally implemented and used within hh-italian-group framework (https://github.com/hh-italian-group/AnalysisTools/).

Configuration is represented as a list of configuration entries. Each configuration entry starts with entry header and ends with an empty line. Entry header has a following format:

```
[TYPE entry_name]
```
*OR*
```
[TYPE entry_name : reference_entry_name]
```
where *TYPE* is entry type and can be omitted for the default entry type, *entry_name* is the entry name that should be unique for a given TYPE, *reference_entry_name* the name of a reference entry of the same type that is used to initialize the entry.

Each configuration entry is a set of properties. Each property is represented by one line in the configuration file and has a following format:
```
property_name: property_value
```
Each entry type defines its own set of properties, related parsing rules and allowed property multiplicity.

Other remarks:
- spaces in entry and property names are not allowed;
- lines that starts with "#" are considered as commentary and are ignored by the configuration parser.

#### StatModel entry properties
Entry type name is *MODEL*. It is the default entry type.

Name                    | Description | Value format | Multiplicity
------------------------|-------------|--------------|-------------
stat_model              | name of the stat model | stat_model_name | &#8804; 1
channels                | list of channels | ch1 ch2 ... | &#8804; 1
categories              | list of categories | cat1 cat2 ... | &#8804; 1
signal_process          | name of the signal process | process_name | &#8804; 1
model_signal_process    | name of the signal process in the format expected for the model dependent interpretation | process_name | &#8804; 1
signal_point_prefix     | prefix before the signal point inside the names of the input shape histograms (e.g. 'M') | prefix | &#8804; 1
signal_points           | list of the signal point (e.g. resonance masses) | p1 p2 ... | &#8804; 1
limit_type              | limit type | model_independent &#124; MSSM | &#8804; 1
th_model_file           | file with the theoretical grid for the model dependent interpretation | file_path | &#8804; 1
blind                   | whatever the final limits should be blind or not | true &#124; false | &#8804; 1
morph                   | whatever the morphing should be applied for input signal shapes (required for model dependent interpretation) | true &#124; false | &#8804; 1
combine_channels        | whatever to combine channels or not | true &#124; false | &#8804; 1
per_channel_limits      | whatever per-channel limits should be produced or not | true &#124; false | &#8804; 1
grid_x                  | range and step definition for the x-axis of the grid points that will be used for model dependent interpretation | min:max:step | &#8804; 1
grid_y                  | range and step definition for the y-axis of the grid points that will be used for model dependent interpretation | min:max:step | &#8804; 1
custom_param            | custom parameter that can be used by the stat model implementation | name value | &#8805; 0

## How to add new statistical model

To add new stat. model one should implement hh_analysis::stat_models::StatModel interface within new class inside HHStatAnalysis/StatModels package and add producer for this class into the producer map in GetProducerFunctions function implemented inside StatModelFactory.cc.
