#!/bin/bash


python scripts/MakeResPlot.py --inputFolder LIMS_LT_Res_500_HMHPC_960_HMMPC_700_LMHPC_960_LMMPC_700_v66/Radion_Node_MASS*/datacards/higgsCombineRadion_Node_MASS --hmFolder LIMS_LT_Res_500_HMHPC_500_HMMPC_000_LMHPC_500_LMMPC_000_v66/Radion_Node_MASS*/datacards/higgsCombineRadion_Node_MASS --lumi 35.9 --log --isAsymptotic --observed --max 3000 --min 0.05 -n ResonantRadionNEW

python scripts/MakeResPlot.py --inputFolder LIMS_LT_Res_500_HMHPC_960_HMMPC_700_LMHPC_960_LMMPC_700_v66/BulkGraviton_Node_MASS*/datacards/higgsCombineBulkGraviton_Node_MASS --hmFolder LIMS_LT_Res_500_HMHPC_500_HMMPC_000_LMHPC_500_LMMPC_000_v66/BulkGraviton_Node_MASS*/datacards/higgsCombineBulkGraviton_Node_MASS --lumi 35.9 --log --isAsymptotic --observed --isGrav --max 3000 --min 0.05 -n ResonantGravitonNEW

#parser.add_argument('--max', dest='max',  default=None, type=float)
#parser.add_argument('--min', dest='min',  default=0.001, type=float)
