// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <string>
#include <utility>
#include <cmath>

//root include files
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

//FLASHgg files
#include "DataFormats/Math/interface/LorentzVector.h"

//Local
#include "HiggsAnalysis/bbggLimits/interface/bbggLTMaker.h"
#include "HiggsAnalysis/bbggLimits/interface/bbgg2DFitter.h"

namespace {
	struct dictionary {
       bbggLTMaker dummy_ltree;
       bbgg2DFitter dummy_fit;
	};
}
